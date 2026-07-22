# ~/~ begin <<docs/src/components/waterdepth.md#src/Components/WaterDepth.jl>>[init]
@compose module WaterDepth
@mixin TimeIntegration, Boxes
using ..Common
using HDF5
using Unitful: ustrip
using ..TimeIntegration: time, time_axis
using CarboKitten.Components.Subsidence
import CarboKitten.Components.Subsidence: cumulative_subsidence

export water_depth, subsider, initial_topography, subsidence_rate_map
export AbstractSubsidenceModifier, MultiplyRate, AddRate, SetRate, Halve, Double
export cumulative_subsidence

# ~/~ begin <<docs/src/components/waterdepth.md#waterdepth-input>>[init]
@kwdef struct Input <: AbstractInput
    sea_level = t -> 0.0u"m"
    initial_topography = (x, y) -> 0.0u"m"
    # scalar Rate (legacy), Matrix{Rate}, or (x,y)->Rate function
    subsidence_rate = 0.0u"m/Myr"
    subsidence_modifiers::Vector{AbstractSubsidenceModifier} = AbstractSubsidenceModifier[]
end
# ~/~ end

# ~/~ begin <<docs/src/components/waterdepth.md#subsidence-rate-map>>[init]
"""
    subsidence_rate_map(input) -> Matrix{Rate}

Normalise `input.subsidence_rate` to a `Matrix{Rate}` of size
`input.box.grid_size`. Accepts a scalar `Quantity`, a `Matrix{Rate}`, or a
function `(x, y) -> Rate`.
"""
function subsidence_rate_map(input::AbstractInput)
    sr = input.subsidence_rate
    if sr isa AbstractMatrix
        @assert size(sr) == input.box.grid_size
        return sr
    elseif sr isa Function
        x, y = box_axes(input.box)
        return sr.(x, y')
    else
        return fill(sr, input.box.grid_size...)
    end
end
# ~/~ end

# ~/~ begin <<docs/src/components/waterdepth.md#subsider>>[init]
function subsider(input::AbstractInput)
    # Scalar legacy path — bit-identical to main branch
    if input.subsidence_rate isa Quantity && isempty(input.subsidence_modifiers)
        Δσ = input.subsidence_rate * input.time.Δt
        return function (state::AbstractState)
            state.bathymetry .-= Δσ
        end
    end

    # Spatio-temporal path — same pattern, per-cell Δσ
    base     = subsidence_rate_map(input)
    x, y     = box_axes(input.box)
    Δt       = input.time.Δt
    t0       = input.time.t0
    get_time = time(input)
    cum      = cumulative_subsidence(base, input.subsidence_modifiers, x, y, t0)

    return function (state::AbstractState)
        t  = get_time(state)
        Δσ = cum(t + Δt) .- cum(t)
        state.bathymetry .-= Δσ
    end
end
# ~/~ end

@kwdef mutable struct State <: AbstractState
    bathymetry::Matrix{Height}
end

@constructor _initial_state(input)::State[bathymetry] =
    (bathymetry = initial_topography(input),)

function initial_state(input::AbstractInput)
    bathymetry = initial_topography(input)
    return State(step=0, bathymetry=bathymetry)
end

function initial_topography(input::AbstractInput)
    if input.initial_topography isa AbstractMatrix
        @assert size(input.initial_topography) == input.box.grid_size
        return input.initial_topography
    end

    x, y = box_axes(input.box)
    return input.initial_topography.(x, y')
end

function water_depth(input::AbstractInput)
    sea_level = input.sea_level
    get_time = time(input)

    return function (state::AbstractState)
        t = get_time(state)
        return sea_level(t) .- state.bathymetry
    end
end

# ~/~ begin <<docs/src/components/waterdepth.md#waterdepth-write-header>>[init]
function write_header(input::AbstractInput, output::AbstractOutput)
    x, y = box_axes(input.box)
    t    = time_axis(input)
    set_attribute(output, "initial_topography", initial_topography(input) |> in_units_of(u"m"))
    set_attribute(output, "sea_level", input.sea_level.(t) .|> in_units_of(u"m"))

    rate_map   = subsidence_rate_map(input)
    scalar_rep = input.subsidence_rate isa Quantity ?
        input.subsidence_rate : sum(rate_map) / length(rate_map)
    set_attribute(output, "subsidence_rate", scalar_rep |> in_units_of(u"m/Myr"))

    if !(input.subsidence_rate isa Quantity)
        set_attribute(output, "subsidence_rate_map", rate_map .|> in_units_of(u"m/Myr"))
    end

    for (idx, m) in enumerate(input.subsidence_modifiers)
        prefix = "subsidence_modifiers/m$(idx)"
        set_attribute(output, "$(prefix)/kind", string(typeof(m).name.name))
        _write_modifier(output, prefix, m)
    end
end

_ser_axis(::Colon)  = [NaN, NaN]
_ser_axis(r::Tuple) = [ustrip(u"m",   r[1]), ustrip(u"m",   r[2])]
_ser_time(::Colon)  = [NaN, NaN]
_ser_time(r::Tuple) = [ustrip(u"Myr", r[1]), ustrip(u"Myr", r[2])]

function _write_modifier(out, prefix, m::MultiplyRate)
    set_attribute(out, "$(prefix)/factor",  m.factor)
    set_attribute(out, "$(prefix)/x_range", _ser_axis(m.x_range))
    set_attribute(out, "$(prefix)/y_range", _ser_axis(m.y_range))
    set_attribute(out, "$(prefix)/t_range", _ser_time(m.t_range))
end
function _write_modifier(out, prefix, m::AddRate)
    set_attribute(out, "$(prefix)/delta",   ustrip(u"m/Myr", m.delta))
    set_attribute(out, "$(prefix)/x_range", _ser_axis(m.x_range))
    set_attribute(out, "$(prefix)/y_range", _ser_axis(m.y_range))
    set_attribute(out, "$(prefix)/t_range", _ser_time(m.t_range))
end
function _write_modifier(out, prefix, m::SetRate)
    set_attribute(out, "$(prefix)/rate",    ustrip(u"m/Myr", m.rate))
    set_attribute(out, "$(prefix)/x_range", _ser_axis(m.x_range))
    set_attribute(out, "$(prefix)/y_range", _ser_axis(m.y_range))
    set_attribute(out, "$(prefix)/t_range", _ser_time(m.t_range))
end
# ~/~ end

end
# ~/~ end
