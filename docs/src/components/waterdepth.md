# Water Depth

```component-dag
CarboKitten.Components.WaterDepth
```

The `WaterDepth` module computes the water depth, given the bedrock elevation, sea level curve, subsidence rate and current sediment height.

## Input

- `initial_topography(x, y)` (a.k.a. initial depth) should be a function taking two coordinates in units of meters, returning an elevation also in meters.
- `sea_level(t)` should be a function taking a time in millions of years (Myr) returning the eustatic sealevel. This could also be an interpolated table.
- `subsidence_rate` the rate of subsidence, which can be specified in three ways:
  - A **scalar** rate (e.g. `50.0u"m/Myr"`) applied uniformly across the grid. This is the legacy behaviour.
  - A **matrix** of rates matching `box.grid_size`, giving each cell its own subsidence rate.
  - A **function** `(x, y) -> Rate` evaluated on the grid axes at construction time.
    Any rate-dimensioned `Quantity` is accepted (`m/yr`, `m/Myr`, etc.).
- `subsidence_modifiers` an optional list of modifiers that locally alter the subsidence rate inside rectangular `(x, y, t)` boxes. See [Subsidence Modifiers](#Subsidence-Modifiers) below.

The signs of these quantities should be such that the following equation holds:

$$T + E = S + W,$$

saying Tectonic subsidence plus Eustatic sea-level change equals Sedimentation plus change in Water depth.

## Subsidence rate map

When subsidence is non-uniform, the helper function `subsidence_rate_map(input)` normalises the input to a `Matrix{Rate}` of size `input.box.grid_size`, regardless of the original specification (scalar, matrix, or function). This mirrors the existing `initial_topography(input)` normaliser.

For backward compatibility, the runtime takes a **scalar fast path** when the subsidence rate is a plain `Quantity` and no modifiers are present — the formula is bit-identical to the original.

## Subsidence Modifiers

Modifiers allow localised changes to the subsidence rate without altering the base rate map. They are defined in the top-level `CarboKitten.Subsidence` module and re-exported by `WaterDepth`. Each modifier is a subtype of `AbstractSubsidenceModifier` with three spatial/temporal range fields:

- `x_range`: `:` (everywhere), or a `(x_lo, x_hi)` tuple.
- `y_range`: `:` (everywhere), or a `(y_lo, y_hi)` tuple.
- `t_range`: `:` (always), or a `(t_lo, t_hi)` tuple.

Inside the box defined by these ranges, the modifier transforms the effective rate. Outside, it has no effect.

### Available modifiers

| Type               | Effect                                     | Key parameter       |
|--------------------|--------------------------------------------|----------------------|
| `MultiplyRate`     | Multiply the rate by a factor              | `factor::Float64`    |
| `AddRate`          | Add a constant to the rate                 | `delta::Rate`        |
| `SetRate`          | Override the rate to a fixed value         | `rate::Rate`         |
| `Halve`            | Shorthand for `MultiplyRate(0.5, ...)`     | —                    |
| `Double`           | Shorthand for `MultiplyRate(2.0, ...)`     | —                    |

Modifiers compose in declaration order: each operates on the rate produced by the previous one.

### Example

```julia
using CarboKitten.Subsidence: MultiplyRate, AddRate, SetRate, Halve

input = ALCAP.Input(
    ...
    subsidence_rate = 50.0u"m/Myr",
    subsidence_modifiers = [
        # Halve subsidence in the first km during [0, 0.5] Myr
        Halve(x_range=(0.0u"m", 1000.0u"m"), t_range=(0.0u"Myr", 0.5u"Myr")),
        # Add 20 m/Myr to a central patch for the entire run
        AddRate(20.0u"m/Myr", x_range=(3.0u"km", 6.0u"km"), y_range=(2.0u"km", 5.0u"km")),
        # Pin subsidence to zero on the right strip for the last quarter
        SetRate(0.0u"m/Myr", x_range=(12.0u"km", 15.0u"km"), t_range=(0.75u"Myr", 1.0u"Myr")),
    ],
    ...)
```

### Cumulative subsidence integral

Given a base rate map $r(x,y)$ and a set of modifiers $\{m_k\}$, the cumulative subsidence at time $t$ is the piecewise integral

$$S(x, y, t) = \int_{t_0}^{t} r_{\mathrm{eff}}(x, y, \tau) \, \mathrm{d}\tau,$$

where $r_{\mathrm{eff}}$ is the base rate sequentially transformed by every modifier whose $(x, y, t)$-box contains the point. Because modifiers are piecewise-constant in time, the integral is computed analytically by walking event boundaries (modifier start/end times intersected with $[t_0, t]$). The closure returned by `cumulative_subsidence(...)` precomputes spatial masks for each modifier, so it is efficient to evaluate at many time steps.

The function `cumulative_subsidence` is generic, with methods dispatching on:

- `(base_rate_map, modifiers, x_axis, y_axis, t0)` → closure, in `Subsidence` module.
- `(header::Header)` → closure, in `Output.Abstract`.
- `(header::Header, t::Time)` → full-grid `Matrix{Length}`.
- `(header::Header, data::Data)` → array matching `data.sediment_thickness` shape.

The last two are used by the plotting routines (sediment profile, Wheeler diagram, glamour view) so that non-uniform subsidence is honoured exactly in post-hoc analysis — not approximated by a scalar mean.

## Source
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
    # Any rate-dimensioned Quantity is accepted (m/yr, m/Myr, etc.).
    subsidence_rate = 0.0u"m/Myr"
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
    input.subsidence_rate isa Quantity && isempty(input.subsidence_modifiers)

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
    scalar_rep = if input.subsidence_rate isa Quantity
        input.subsidence_rate
    else
        sum(rate_map) / length(rate_map)
    end
    set_attribute(output, "subsidence_rate", scalar_rep |> in_units_of(u"m/Myr"))

    # Per-cell rate map (only when non-uniform). The AbstractArray dispatch of
    # `set_attribute` writes a dataset, so this lands under `input/subsidence_rate_map`.
    if !(input.subsidence_rate isa Quantity)
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

## Subsidence module

The modifier types and cumulative-subsidence algorithm are defined in the [Subsidence](subsidence.md) module. See that page for the full source and API.
