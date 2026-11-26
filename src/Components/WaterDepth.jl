# ~/~ begin <<docs/src/components/waterdepth.md#src/Components/WaterDepth.jl>>[init]
@component module WaterDepth
@mixin TimeIntegration, Boxes
using ..Common
import ..Common: write_header
using HDF5
using ..TimeIntegration: time, time_axis

export water_depth, subsider, initial_topography, cumulative_subsidence
export TimeVariableSubsidenceRate, SpaceTimeVariableSubsidenceRate

@kwdef struct Input <: InputPart
    sea_level = t -> 0.0u"m"
    initial_topography = (x, y) -> 0.0u"m"
    subsidence_rate = 0.0u"m/Myr"
end

@kwdef mutable struct State <: StatePart
    bathymetry::Matrix{Height}
end

Common.initial_state(::Input, input::AbstractInput) =
    State(bathymetry=initial_topography(input))

initial_state(input::AbstractInput) = Common.initial_state(input)

function initial_topography(input::AbstractInput)
    if input.initial_topography isa AbstractMatrix
        @assert size(input.initial_topography) == input.box.grid_size
        return input.initial_topography
    end

    x, y = box_axes(input.box)
    return input.initial_topography.(x, y')
end

struct TimeVariableSubsidenceRate
    f
end

struct SpaceTimeVariableSubsidenceRate
    f
end

subsider(input::AbstractInput) = subsider(input, input.subsidence_rate)

function subsider(input::AbstractInput, s::SpaceTimeVariableSubsidenceRate)
    get_time = time(input)
    x, y = box_axes(input.box)
    dt = input.time.Δt

    return function (state::AbstractState)
        t = get_time(state)
        state.bathymetry .-= s.f.(x, y', t) .* dt
    end
end

function subsider(input::AbstractInput, s::TimeVariableSubsidenceRate)
    get_time = time(input)
    dt = input.time.Δt

    return function (state::AbstractState)
        t = get_time(state)
        Δσ = s.f(t) * dt
        state.bathymetry .-= Δσ
    end
end

function subsider(input::AbstractInput, s::Quantity)
    Δσ = input.subsidence_rate * input.time.Δt

    return function (state::AbstractState)
        state.bathymetry .-= Δσ
    end
end

"""
    cumulative_subsidence(input)

Total amount of subsidence over the entire run, for each grid cell.
"""
function cumulative_subsidence(input::AbstractInput)
    state = record_type(:State, WaterDepth)(step=0, bathymetry=zeros(Amount, input.box.grid_size...))
    subside! = subsider(input)
    for _ in 1:input.time.steps
        subside!(state)
        state.step += 1
    end
    return .-state.bathymetry
end

"""
    water_depth(input) -> f

Returns a function `f(state)` that computes the water depth for the given state.
The result is written into a buffer owned by `f`, so the returned array is
overwritten on the next call.
"""
function water_depth(input::AbstractInput)
    sea_level = input.sea_level
    get_time = time(input)
    out = Matrix{Height}(undef, input.box.grid_size...)

    return function (state::AbstractState)
        t = get_time(state)
        out .= sea_level(t) .- state.bathymetry
        return out
    end
end

function write_header(::Input, input::AbstractInput, output::AbstractOutput)
    x, y = box_axes(input.box)
    t = time_axis(input)
    set_attribute(output, "initial_topography", initial_topography(input) |> in_units_of(u"m"))
    set_attribute(output, "sea_level", input.sea_level.(t) .|> in_units_of(u"m"))
    set_attribute(output, "total_subsidence", cumulative_subsidence(input) |> in_units_of(u"m"))
    if input.subsidence_rate isa Quantity
        set_attribute(output, "subsidence_rate", input.subsidence_rate |> in_units_of(u"m/Myr"))
    end
end

end
# ~/~ end
