# Water Depth

```component-dag
CarboKitten.Components.WaterDepth
```

The `WaterDepth` module computes the water depth, given the bedrock elevation, sea level curve, subsidence rate and current sediment height.

## Input

- `initial_topography(x, y)` (a.k.a. initial depth) should be a function taking two coordinates in units of meters, returning an elevation also in meters.
- `sea_level(t)` should be a function taking a time in millions of years (Myr) returning the eustatic sealevel. This could also be an interpolated table.
- `subsidence_rate` a constant rate of subsidence in m/Myr.

The signs of these quantities should be such that the following equation holds:

$$T + E = S + W,$$

saying Tectonic subsidence plus Eustatic sea-level change equals Sedimentation plus change in Water depth.

``` {.julia file=src/Components/WaterDepth.jl}
@compose module WaterDepth
@mixin TimeIntegration, Boxes
using ..Common
using HDF5
using Unitful
using Unitful: ustrip
using ..TimeIntegration: time, time_axis
using CarboKitten.Subsidence
import CarboKitten.Subsidence: cumulative_subsidence

export water_depth, subsidence_rate_map
export AbstractSubsidenceModifier, MultiplyRate, AddRate, SetRate, Halve, Double, apply_rate
export cumulative_subsidence

@kwdef struct Input <: AbstractInput
    sea_level = t -> 0.0u"m"
    initial_topography = (x, y) -> 0.0u"m"
    # Either a scalar `Rate` (legacy) or a `Matrix{Rate}` matching `box.grid_size`,
    # or a `(x, y) -> Rate` function evaluated on the grid at construction time.
    subsidence_rate::Union{Rate,AbstractMatrix{<:Rate},Function} = 0.0u"m/Myr"
    subsidence_modifiers::Vector{AbstractSubsidenceModifier} = AbstractSubsidenceModifier[]
    end

@kwdef mutable struct State <: AbstractState
    sediment_height::Matrix{Height}
end

function initial_state(input::AbstractInput)
    return State(step=0, sediment_height=zeros(Height, input.box.grid_size...))
end

function initial_topography(input::AbstractInput)
    if input.initial_topography isa AbstractMatrix
        @assert size(input.initial_topography) == input.box.grid_size
        return input.initial_topography
    end

    x, y = box_axes(input.box)
    return input.initial_topography.(x, y')
end

"""
    subsidence_rate_map(input) -> Matrix{Rate}

Return the base subsidence rate as a `Matrix{Rate}` of size `input.box.grid_size`,
regardless of whether `input.subsidence_rate` was given as a scalar, a matrix, or
a function. Mirrors the `initial_topography(input)` normalizer.
"""
function subsidence_rate_map(input::AbstractInput)
    sr = input.subsidence_rate
    if sr isa AbstractMatrix
        @assert size(sr) == input.box.grid_size "subsidence_rate matrix size $(size(sr)) does not match box.grid_size $(input.box.grid_size)"
        return sr
    elseif sr isa Function
        x, y = box_axes(input.box)
        return sr.(x, y')
    else
        return fill(sr, input.box.grid_size...)
    end
end

# Scalar shortcut: if subsidence is uniform and there are no modifiers, the
# legacy formula is bit-identical and avoids a temporary matrix per timestep.
_is_scalar_subsidence(input::AbstractInput) =
    input.subsidence_rate isa Rate && isempty(input.subsidence_modifiers)

"""
    cumulative_subsidence(input::AbstractInput) -> (t -> Matrix{Length})

`AbstractInput` method of `CarboKitten.Subsidence.cumulative_subsidence`.
Returns a closure that, given absolute time `t`, returns the per-cell
cumulative subsidence amount integrated from `input.time.t0`.
"""
function cumulative_subsidence(input::AbstractInput)
    base = subsidence_rate_map(input)
    x, y = box_axes(input.box)
    return cumulative_subsidence(base, input.subsidence_modifiers, x, y, input.time.t0)
end

function water_depth(input::AbstractInput)
    x, y = box_axes(input.box)
    eta0 = initial_topography(input)
    sea_level = input.sea_level
    t0 = input.time.t0
    get_time = time(input)

    if _is_scalar_subsidence(input)
        # Legacy bit-identical path
        sr = input.subsidence_rate
        return function (state::AbstractState)
            t = get_time(state)
            return sea_level(t) .- eta0 .+
                   (sr * (t - t0)) .- state.sediment_height
        end
    else
        cum = cumulative_subsidence(input)
        return function (state::AbstractState)
            t = get_time(state)
            return sea_level(t) .- eta0 .+
                   cum(t) .- state.sediment_height
        end
    end
end

function write_header(input::AbstractInput, output::AbstractOutput)
    x, y = box_axes(input.box)
    t = time_axis(input)
    set_attribute(output, "initial_topography", initial_topography(input) |> in_units_of(u"m"))
    set_attribute(output, "sea_level", input.sea_level.(t) .|> in_units_of(u"m"))

    # Backward-compatible scalar summary so old consumers reading
    # `header.subsidence_rate` keep working.
    rate_map = subsidence_rate_map(input)
    scalar_rep = if input.subsidence_rate isa Rate
        input.subsidence_rate
    else
        sum(rate_map) / length(rate_map)
    end
    set_attribute(output, "subsidence_rate", scalar_rep |> in_units_of(u"m/Myr"))

    # Per-cell rate map (only when non-uniform). The AbstractArray dispatch of
    # `set_attribute` writes a dataset, so this lands under `input/subsidence_rate_map`.
    if !(input.subsidence_rate isa Rate)
        set_attribute(output, "subsidence_rate_map",
                      rate_map .|> in_units_of(u"m/Myr"))
    end

    # Modifier descriptors, one group per modifier.
    for (idx, m) in enumerate(input.subsidence_modifiers)
        prefix = "subsidence_modifiers/m$(idx)"
        set_attribute(output, "$(prefix)/kind", string(typeof(m).name.name))
        _serialize_modifier(output, prefix, m)
    end
end

_serialize_axis_spec(::Colon)            = [NaN, NaN]
_serialize_axis_spec(r::Tuple)           = [ustrip(u"m",   r[1]), ustrip(u"m",   r[2])]
_serialize_axis_spec(r::AbstractRange)   = [ustrip(u"m",   first(r)), ustrip(u"m",   last(r))]
_serialize_time_spec(::Colon)            = [NaN, NaN]
_serialize_time_spec(r::Tuple)           = [ustrip(u"Myr", r[1]), ustrip(u"Myr", r[2])]
_serialize_time_spec(r::AbstractRange)   = [ustrip(u"Myr", first(r)), ustrip(u"Myr", last(r))]

function _serialize_modifier(output, prefix, m::MultiplyRate)
    set_attribute(output, "$(prefix)/factor",  m.factor)
    set_attribute(output, "$(prefix)/x_range", _serialize_axis_spec(m.x_range))
    set_attribute(output, "$(prefix)/y_range", _serialize_axis_spec(m.y_range))
    set_attribute(output, "$(prefix)/t_range", _serialize_time_spec(m.t_range))
end
function _serialize_modifier(output, prefix, m::AddRate)
    set_attribute(output, "$(prefix)/delta",   ustrip(u"m/Myr", m.delta))
    set_attribute(output, "$(prefix)/x_range", _serialize_axis_spec(m.x_range))
    set_attribute(output, "$(prefix)/y_range", _serialize_axis_spec(m.y_range))
    set_attribute(output, "$(prefix)/t_range", _serialize_time_spec(m.t_range))
end
function _serialize_modifier(output, prefix, m::SetRate)
    set_attribute(output, "$(prefix)/rate",    ustrip(u"m/Myr", m.rate))
    set_attribute(output, "$(prefix)/x_range", _serialize_axis_spec(m.x_range))
    set_attribute(output, "$(prefix)/y_range", _serialize_axis_spec(m.y_range))
    set_attribute(output, "$(prefix)/t_range", _serialize_time_spec(m.t_range))
end

end
```
