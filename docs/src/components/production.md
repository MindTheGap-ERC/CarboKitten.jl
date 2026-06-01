# Production

```component-dag
CarboKitten.Components.Production
```

The `Production` module specifies the production rate following the model by Bosscher & Schlager 1992 [Bosscher1992](@cite).
The growth rate is given as

$$g(w) = g_m \tanh\left({{I_0 e^{-kw}} \over {I_k}}\right).$$

This can be understood as a smooth transition between the maximum growth rate under saturated conditions, and exponential decay due to light intensity dropping with greater water depth.

``` {.julia #component-production-rate}
function production_rate(insolation, facies, water_depth)
    gₘ = facies.maximum_growth_rate
    I = insolation / facies.saturation_intensity
    x = water_depth * facies.extinction_coefficient
    return x > 0.0 ? gₘ * tanh(I * exp(-x)) : zero(typeof(gₘ))
end

benthic_production(i, f, w) = production_rate(i, f, w)

function capped_production(f, insolation, water_depth, dt, factor::Real=1.0)
    clip_positive(x::T) where {T} = max(x, zero(T))
    p = clip_positive(f(insolation, water_depth)) * factor
    return min(max(0.0u"m", water_depth), p * dt)
end
```

The `capped_production` function accepts an optional `factor` (default `1.0`) that scales the production rate before capping. This is used by the [time-window modifiers](#Time-window-production-modifiers) to vary production over time without changing the base curve. Existing call sites with four arguments are unaffected.

From just this equation we can define a uniform production process. This requires that we have a `Facies` that defines the `maximum_growth_rate`, `extinction_coefficient` and `saturation_intensity`.

The `insolation` input may be given as a scalar quantity, say `400u"W/m^2"`, or as a function of time.

## Pelagic Production

A facies can be specified as being pelagic, meaning that production is not governed at the sea floor (benthic zone), rather the entire water column above contributes to the production. The local production rate still follows the same function as before (well motivated by exponential decay of available radiation), but now we need to integrate over the water depth:

$$g_p(w) = \int_0^{w} g_m \tanh\left({{I_0 e^{-kw}} \over {I_k}}\right) \textrm{d}w.$$

This integral has no analytic solution, so we'll use a numeric integrator to evaluate the complete production curve, and then linearly interpolate the generated table to obtain production rates.

Pelagic facies do not participate in the CA.

``` {.julia #pelagic-production}
function pelagic_production(insolation, facies, water_depth)
    return quadgk(w -> production_rate(insolation, facies, w), 0.0u"m", water_depth)[1]
end
```

``` {.julia file=examples/production/pelagic.jl}
module PelagicProductionPlot

using CarboKitten
using CarboKitten.Production
using CairoMakie

@kwdef struct Input <: CarboKitten.AbstractInput
    insolation = 400.0u"W/m^2"
end

function main()
    water_depth = (0.01:0.1:50.0)u"m"
    fig = Figure(size=(600,600))
    input = Input()

    ax = Axis(fig[1, 1], yreversed=true, ylabel="depth [m]", xlabel="production [m/Myr]")
    for (k, prod) in pairs(Production.EXAMPLE)
        f = production_profile(input, prod)
        p = water_depth .|> (w -> f(input.insolation, w))
        lines!(ax, p |> in_units_of(u"m/Myr"),
            water_depth |> in_units_of(u"m"), label = string(k))
    end
    fig[1, 2] = Legend(fig, ax, "Facies")
    fig
end

end
```

## Interpolated Production

Not all production curves follow the Bosscher & Schlager model. In some cases it is more natural to specify the curve directly, as a set of depth knots and multipliers applied to a peak rate. The `InterpolatedProduction` type supports this:

$$g(w) = g_{\max} \cdot f(w),$$

where $f(w)$ is a piecewise-linear function defined by `(depth_knots, multipliers)` pairs, with **flat extrapolation** outside the knot range (i.e. the rate at $w < \min(\text{knots})$ equals $g_{\max} \cdot m_{\min}$, and similarly at the upper end).

This curve is **independent of insolation** — the per-depth shape is fixed by the user. To vary the overall scale over time, pair it with the [time-window modifiers](#Time-window-production-modifiers) described below.

### Example

```julia
using CarboKitten.Production: InterpolatedProduction

# Shallow reef builder: peaks at 5–15 m, dies off by 50 m
InterpolatedProduction(
    maximum_production = 500.0u"m/Myr",
    depth_knots  = [0.0u"m", 5.0u"m", 15.0u"m", 30.0u"m", 50.0u"m"],
    multipliers  = [0.0,     1.0,     1.0,      0.4,      0.0])
```

At 10 m depth the rate is `500 × 1.0 = 500 m/Myr`; at 30 m it is `500 × 0.4 = 200 m/Myr`.

The `depth_knots` need not be sorted — they are sorted internally when building the interpolator. `length(depth_knots)` must equal `length(multipliers)` and be ≥ 2.

### Mixing curve types

Different facies in the same run can use different production types. For example:

```julia
facies = [
    ALCAP.Facies(production = BenthicProduction(
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2")),
    ALCAP.Facies(production = InterpolatedProduction(
        maximum_production=350.0u"m/Myr",
        depth_knots=[0.0u"m", 10.0u"m", 40.0u"m", 60.0u"m"],
        multipliers=[0.0, 0.5, 1.0, 0.0])),
]
```

The `production_curve` visualisation handles all three types (benthic, pelagic, interpolated) both from an in-memory input and from an HDF5 file.

## Time-window production modifiers

Production modifiers scale the production rate of selected facies during a specified time window. They are defined alongside the production types in `src/Production.jl` and composed multiplicatively in declaration order.

### `MultiplyProduction`

The only modifier type currently provided. It multiplies the production rate by `factor` for `facies` during `t_range`:

```julia
using CarboKitten.Production: MultiplyProduction

production_modifiers = [
    # Halve all facies during [0, 0.3] Myr
    MultiplyProduction(0.5; t_range=(0.0u"Myr", 0.3u"Myr")),
    # Double facies 1 during [0.5, 0.8] Myr
    MultiplyProduction(2.0; t_range=(0.5u"Myr", 0.8u"Myr"), facies=1),
]
```

Parameters:

- `factor::Float64` — multiplicative scaling factor.
- `t_range` — `:` (always active) or a `(t_lo, t_hi)` tuple in time units.
- `facies` — `:` (all facies), a single `Int` (1-based index), or a `Vector{Int}`.

The net factor at time $t$ for facies $k$ is:

$$f_{\mathrm{eff}}(t, k) = \prod_{m \in \mathrm{active}(t, k)} m.\mathrm{factor}$$

### How it works

At each time step, the runtime production functions (`uniform_production` in `Components.Production` and `production` in `CAProduction`) compute a per-facies factor vector by calling `production_factor(modifiers, t, facies_index)`. This factor is passed as the fifth argument to `capped_production`, scaling the rate before the cap is applied. When no modifiers are present, the factor is `1.0` — the legacy code path is bit-identical.

## Production Component

### Production Properties
Because the parameters for benthic and pelagic production have different units, we need different types to store them.

``` {.julia #production-profile}
"""
    production_profile(input::AbstractInput, p)

Given an input production configuration, returns a function of
insolation (in W/m^2) and water depth (in m) to production (in m/Myr).

The default implementation is the identity function; it assumes
that its parameter is a callable object with the correct behaviour.

## Example
If you want to implement your own production profile, define a new type,

    struct MyProduction <: AbstractProduction
        ...
    end

and implement this method for that type,

    import CarboKitten.Components.Production: production_profile

    production_profile(input::AbstractInput, p::MyProduction) =
        function (insolation, water_depth)
            ...
        end
"""
production_profile(input::AbstractInput, p) = p

"""
    is_benthic(obj)

Predicate to determine if a facies or production spec is benthic.
Defaults to `false`.
"""
is_benthic(p) = false

"""
    is_pelagic(obj)

Predicate to determine if a facies or production spec is pelagic.
Defaults to `false`.
"""
is_pelagic(p) = false

abstract type AbstractProduction end

struct NoProduction <: AbstractProduction
end

production_profile(::AbstractInput, ::NoProduction) = (_, _) -> 0.0u"m/Myr"

@kwdef struct BenthicProduction <: AbstractProduction
    maximum_growth_rate::typeof(1.0u"m/Myr") = 0.0u"m/Myr"
    extinction_coefficient::typeof(1.0u"m^-1") = 0.0u"m^-1"
    saturation_intensity::typeof(1.0u"W/m^2") = 1.0u"W/m^2"
end

is_benthic(::BenthicProduction) = true
is_pelagic(::BenthicProduction) = false
production_profile(input::AbstractInput, p::BenthicProduction) = (I, w) -> benthic_production(I, p, w)

@kwdef struct PelagicProduction <: AbstractProduction
    maximum_growth_rate::typeof(1.0u"1/Myr") = 0.0u"1/Myr"
    extinction_coefficient::typeof(1.0u"m^-1") = 0.0u"m^-1"
    saturation_intensity::typeof(1.0u"W/m^2") = 1.0u"W/m^2"
    maximum_production_depth::typeof(1.0u"m") = 200.0u"m"
    table_size::Tuple{Int, Int} = (1000, 1000)
end

is_benthic(::PelagicProduction) = false
is_pelagic(::PelagicProduction) = true
production_profile(input::AbstractInput, p::PelagicProduction) = pelagic_production_lookup(input, p)

```

### Input parameters

In the case of pelagic production, the exact production rate can be quite expensive to compute. Instead, we create look-up tables for each facies. These tables are also easy to store in HDF5. The `maximum_production_depth` parameter controls down to which depth the table is precomputed. At larger depths, we use linear extrapolation (negative values will be clipped to 0). The `production_table_size` parameter determines the resolution of the look-up table.

``` {.julia #production-input}
@kwdef struct Input <: AbstractInput
    insolation
    production_modifiers = []
end

@kwdef struct Facies <: AbstractFacies
    production = NoProduction()
end

is_benthic(facies::AbstractFacies) = is_benthic(facies.production)
is_pelagic(facies::AbstractFacies) = is_pelagic(facies.production)
```

The `production_modifiers` field defaults to an empty list, so existing inputs that don't specify it keep their original behaviour. It is left untyped because the `@compose` macro processes `@kwdef` struct fields at expansion time, before `using` statements have executed — typed references to `AbstractProductionModifier` would fail at that stage.

### Insolation

Insolation can be passed as a constant, function or array (which should have the same length as number of time steps in the run).

``` {.julia #production-insolation}
function insolation(input::AbstractInput)
    insolation = input.insolation
    tprop = input.time
    if insolation isa Quantity
        return s -> insolation
    end
    if insolation isa AbstractVector
        @info "Reading insolation from a table"
        return s -> insolation[s.step+1]
    end
    @info "Reading insolation from a function"
    function (s::AbstractState)
        t = time(tprop, s)
        return insolation(t)
    end
end
```

### Pelagic Lookup tables

We use `linear_interpolation` from `Interpolations` to compute production profiles from look-up tables. If insolation is constant we can use a one-dimensional interpolation. In case of variable insolation, we first find the extrema of the input insolation curve, and then compute a two-dimensional table of values.

``` {.julia #production-lookup}
function pelagic_production_lookup(input::AbstractInput, prod::PelagicProduction)
    insolation_param = input.insolation
    depth_grid = LinRange(0.0, prod.maximum_production_depth |> in_units_of(u"m"), prod.table_size[2])

    if insolation_param isa Quantity
        production_values = [pelagic_production(insolation_param, prod, w * u"m") |> in_units_of(u"m/Myr") for w in depth_grid]
        interpolated = linear_interpolation(depth_grid, production_values, extrapolation_bc = Line())
        return (I0, w) -> interpolated(w |> in_units_of(u"m")) * u"m/Myr"
    end

    insolation_vec = if insolation_param isa AbstractVector
        insolation_param
    else
        t = time_axis(input)
        insolation_param.(t)
    end
    I_min, I_max = extrema(insolation_vec)

    insolation_grid = LinRange(
        I_min |> in_units_of(u"W/m^2"),
        I_max |> in_units_of(u"W/m^2"), prod.table_size[1])
    production_values = [
        pelagic_production(I, prod, w * u"m") |> in_units_of(u"m/Myr")
        for I in insolation_grid, w in depth_grid
    ]
    interpolated = linear_interpolation((insolation_grid, depth_grid), production_values, extrapolation_bc = Line())
    return (I0, w) -> interpolated(I0 |> in_units_of(u"W/m^2"), w |> in_units_of(u"m")) * u"m/Myr"
end
```

### HDF5 serialization

Production types and modifiers are serialized under the `input/` group:

- **Benthic/Pelagic facies**: `facies$(i)/type` attribute (`"benthic"` or `"pelagic"`), plus `maximum_growth_rate`, `extinction_coefficient`, `saturation_intensity` as attributes.
- **Interpolated facies**: `facies$(i)/type` attribute = `"interpolated"`, plus:
  - `maximum_production` — scalar attribute in `m/Myr`.
  - `depth_knots` — dataset (array) in meters.
  - `multipliers` — dataset (array) of `Float64`.

  Note that `depth_knots` and `multipliers` are written as **datasets** (via the `AbstractArray` dispatch of `set_attribute`), not as scalar attributes. The HDF5 reader in `ProductionCurve.jl` reads them from the group directly, not from `HDF5.attributes`.

- **Production modifiers**: one group per modifier under `input/production_modifiers/m1/`, `.../m2/`, etc., each containing:
  - `kind` — string attribute (`"MultiplyProduction"`).
  - `factor` — scalar `Float64` attribute.
  - `t_range` — 2-element array in Myr; `[NaN, NaN]` encodes `:`.
  - `facies` — integer array; `[-1]` encodes `:` (all facies).

Old HDF5 files without these fields are read normally; missing facies type attributes default to the existing benthic/pelagic logic.

### Boilerplate

We put the basic production equations in a separate module `CarboKitten.Production`. This will contain the dynamic API (`production_profile`) and the implementation of `BenthicProduction` and `PelagicProduction`. We isolate these from the component `CarboKitten.Components.Production` implementation to prevent these production objects from being replicated in derived components, leading to problems with dispatch on `production_profile`.

This design makes sense, since we isolate production equations from the larger lego-brick design of CarboKitten.

``` {.julia file=src/Production.jl}
module Production

using Unitful
using QuadGK
using Interpolations

import ..CarboKitten: AbstractInput, time_axis
using ..Utility: in_units_of

export production_profile

<<component-production-rate>>
<<pelagic-production>>
<<production-profile>>
<<production-lookup>>

# =============================================================================
# Interpolated (knot-based) production curve
# =============================================================================

"""
    InterpolatedProduction(; maximum_production, depth_knots, multipliers)

A depth-only production curve defined by a peak rate and a piecewise-linear
shape over `(depth, multiplier)` knots:

    rate(w) = maximum_production × interpolate(depth_knots, multipliers; w)

Linear interpolation between knots, flat extrapolation outside (closest knot
value). Independent of insolation.

`depth_knots` need not be sorted; they are sorted internally. `length(depth_knots)`
must equal `length(multipliers)` and be ≥ 2.

# Example

```julia
InterpolatedProduction(
    maximum_production = 500.0u"m/Myr",
    depth_knots        = [0.0u"m", 5.0u"m", 20.0u"m", 50.0u"m"],
    multipliers        = [0.0,     1.0,     0.6,      0.0])
```
"""
@kwdef struct InterpolatedProduction <: AbstractProduction
    maximum_production::typeof(1.0u"m/Myr") = 0.0u"m/Myr"
    depth_knots::Vector{typeof(1.0u"m")}    = typeof(1.0u"m")[]
    multipliers::Vector{Float64}            = Float64[]
end

is_benthic(::InterpolatedProduction) = true
is_pelagic(::InterpolatedProduction) = false

function production_profile(::AbstractInput, p::InterpolatedProduction)
    @assert length(p.depth_knots) == length(p.multipliers) "InterpolatedProduction: depth_knots and multipliers must have equal length"
    @assert length(p.depth_knots) >= 2 "InterpolatedProduction: need at least 2 knots"
    depths_m = [d |> in_units_of(u"m") for d in p.depth_knots]
    order = sortperm(depths_m)
    itp = linear_interpolation(depths_m[order], p.multipliers[order], extrapolation_bc=Flat())
    max_rate = p.maximum_production
    return (_, w) -> max_rate * itp(w |> in_units_of(u"m"))
end

# =============================================================================
# Time-window production modifiers
# =============================================================================

const _ProdTime      = typeof(1.0u"Myr")
const _ProdTimeSpec  = Union{Colon, Tuple{_ProdTime,_ProdTime}}
const _FaciesSpec    = Union{Colon, Int, AbstractVector{Int}}

abstract type AbstractProductionModifier end

"""
    MultiplyProduction(factor; t_range=:, facies=:)

Multiply the production rate by `factor` during `t_range` for `facies`.
`t_range` is `:` (always) or a `(t_lo, t_hi)` tuple. `facies` is `:`
(all), a single `Int`, or a `Vector{Int}`.
"""
@kwdef struct MultiplyProduction <: AbstractProductionModifier
    factor::Float64
    t_range::_ProdTimeSpec = (:)
    facies::_FaciesSpec    = (:)
end
MultiplyProduction(factor::Real; kwargs...) = MultiplyProduction(; factor=Float64(factor), kwargs...)

"""Net multiplicative factor from all active modifiers at time `t` for facies index `fi`."""
function production_factor(mods, t, fi::Int)
    f = 1.0
    for m in mods
        active_t = m.t_range isa Colon || (m.t_range[1] <= t <= m.t_range[2])
        active_f = m.facies isa Colon || (m.facies isa Int ? m.facies == fi : fi in m.facies)
        if active_t && active_f
            f *= m.factor
        end
    end
    return f
end

const EXAMPLE = Dict(
    :euphotic => BenthicProduction(
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2"),
    :oligophotic => BenthicProduction(
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2"),
    :aphotic => BenthicProduction(
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2"),
    :pelagic => PelagicProduction(
        maximum_growth_rate=7.0u"1/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2"),
    :interpolated => InterpolatedProduction(
        maximum_production=500u"m/Myr",
        depth_knots=[0.0u"m", 5.0u"m", 15.0u"m", 30.0u"m", 50.0u"m"],
        multipliers=[0.0, 1.0, 1.0, 0.4, 0.0])
)

end
```

``` {.julia file=src/Components/Production.jl}
@compose module Production
@mixin TimeIntegration, WaterDepth, FaciesBase
using ..Common
using ..WaterDepth: water_depth
using ..TimeIntegration: time, write_times
using ...Production: NoProduction, InterpolatedProduction,
    AbstractProductionModifier, MultiplyProduction, production_factor
import ...Production: production_profile, is_benthic, is_pelagic, capped_production

using HDF5
using QuadGK
using Interpolations
using Logging

export uniform_production
export AbstractProductionModifier, MultiplyProduction, InterpolatedProduction, NoProduction

<<production-input>>
<<production-insolation>>

function write_header(input::AbstractInput, output::AbstractOutput)
    if input.insolation isa Quantity
        set_attribute(output, "insolation",
            fill(input.insolation |> in_units_of(u"W/m^2"), input.time.steps))
    elseif input.insolation isa AbstractVector
        set_attribute(output, "insolation",
            input.insolation |> in_units_of(u"W/m^2"))
    else
        t = write_times(input)[1:end-1]
        set_attribute(output, "insolation",
            input.insolation.(t) |> in_units_of(u"W/m^2"))
    end

    for (i, f) in enumerate(input.facies)
        p = f.production
        if p isa PelagicProduction
            set_attribute(output, "facies$(i)/type", "pelagic")
            set_attribute(output, "facies$(i)/maximum_growth_rate", p.maximum_growth_rate |> in_units_of(u"1/Myr"))
            set_attribute(output, "facies$(i)/extinction_coefficient", p.extinction_coefficient |> in_units_of(u"m^-1"))
            set_attribute(output, "facies$(i)/saturation_intensity", p.saturation_intensity |> in_units_of(u"W/m^2"))
        elseif p isa InterpolatedProduction
            set_attribute(output, "facies$(i)/type", "interpolated")
            set_attribute(output, "facies$(i)/maximum_production", p.maximum_production |> in_units_of(u"m/Myr"))
            set_attribute(output, "facies$(i)/depth_knots", p.depth_knots .|> in_units_of(u"m"))
            set_attribute(output, "facies$(i)/multipliers", collect(Float64, p.multipliers))
        elseif is_benthic(p)
            set_attribute(output, "facies$(i)/type", "benthic")
            set_attribute(output, "facies$(i)/maximum_growth_rate", p.maximum_growth_rate |> in_units_of(u"m/Myr"))
            set_attribute(output, "facies$(i)/extinction_coefficient", p.extinction_coefficient |> in_units_of(u"m^-1"))
            set_attribute(output, "facies$(i)/saturation_intensity", p.saturation_intensity |> in_units_of(u"W/m^2"))
        end
    end

    for (idx, m) in enumerate(input.production_modifiers)
        prefix = "production_modifiers/m$(idx)"
        set_attribute(output, "$(prefix)/kind", "MultiplyProduction")
        set_attribute(output, "$(prefix)/factor", m.factor)
        set_attribute(output, "$(prefix)/t_range",
            m.t_range isa Colon ? [NaN, NaN] : [ustrip(u"Myr", m.t_range[1]), ustrip(u"Myr", m.t_range[2])])
        set_attribute(output, "$(prefix)/facies",
            m.facies isa Colon ? [-1] : m.facies isa Int ? [m.facies] : collect(Int, m.facies))
    end
end

function uniform_production(input::AbstractInput)
    w = water_depth(input)
    na = [CartesianIndex()]
    insolation_func = insolation(input)
    facies = input.facies
    dt = input.time.Δt
    production_rates = [production_profile(input, f.production) for f in facies]

    modifiers = hasfield(typeof(input), :production_modifiers) ?
        input.production_modifiers : AbstractProductionModifier[]

    if isempty(modifiers)
        p(state::AbstractState, wd::AbstractMatrix) =
            capped_production.(production_rates[:, na, na], insolation_func(state), wd[na, :, :], dt)
        p(state::AbstractState) = p(state, w(state))
        return p
    end

    n_f = length(facies)
    get_time = time(input)
    function p_mod(state::AbstractState, wd::AbstractMatrix)
        t = get_time(state)
        factors = Float64[production_factor(modifiers, t, i) for i in 1:n_f]
        return capped_production.(production_rates[:, na, na], insolation_func(state),
                                  wd[na, :, :], dt, factors[:, na, na])
    end
    p_mod(state::AbstractState) = p_mod(state, w(state))
    return p_mod
end

end
```

Note the ordering of `isa` checks in `write_header`: `InterpolatedProduction` must be checked **before** `is_benthic(p)` because `is_benthic(::InterpolatedProduction) = true` — without the `isa` guard, the benthic branch would try to read `maximum_growth_rate`, which `InterpolatedProduction` does not have.

## CA Production

```component-dag
CarboKitten.Components.CAProduction
```

The `CAProduction` component gives production that depends on the provided CA. The production modifier mechanism is threaded through identically: a per-facies factor tuple is computed once per step and passed to `capped_production`.

``` {.julia file=src/Components/CAProduction.jl}
@compose module CAProduction
    @mixin TimeIntegration, CellularAutomaton, Production
    using ..Common
    using ..Production: insolation
    using ..TimeIntegration: time
    using ..WaterDepth: water_depth
    using ...Production: production_profile, capped_production,
        AbstractProductionModifier, production_factor
    using Logging

    function production(input::AbstractInput)
        w = water_depth(input)
        na = [CartesianIndex()]
        output_ = Array{Amount, 3}(undef, n_facies(input), input.box.grid_size...)

        w = water_depth(input)
        s = insolation(input)
        n_f = n_facies(input)
        facies = input.facies
        active_facies = [i for (i, f) in pairs(facies) if f.active]
        global_facies = [i for (i, f) in pairs(facies) if !f.active]
        dt = input.time.Δt
        # Having this a Tuple should make things type stable?
        production_specs = ((production_profile(input, f.production) for f in facies)...,)

        modifiers = hasfield(typeof(input), :production_modifiers) ?
            input.production_modifiers : AbstractProductionModifier[]
        has_mods = !isempty(modifiers)
        get_time = has_mods ? time(input) : nothing

        function p(state::AbstractState, wd::AbstractMatrix)::Array{Amount,3}
            output::Array{Amount, 3} = output_
            insolation::typeof(1.0u"W/m^2") = s(state)
            factors = if has_mods
                t = get_time(state)
                ntuple(i -> production_factor(modifiers, t, i), n_f)
            else
                ntuple(_ -> 1.0, n_f)
            end
            for i in eachindex(IndexCartesian(), wd)
                for f in eachindex(facies)
                    if facies[f].active
                        output[f, i[1], i[2]] = f != state.ca[i] ? 0.0u"m" :
                            capped_production(production_specs[f], insolation, wd[i], dt, factors[f])
                    else
                        output[f, i[1], i[2]] =
                            capped_production(production_specs[f], insolation, wd[i], dt, factors[f])
                    end
                end
            end
            return output
        end

        @inline p(state::AbstractState) = p(state, w(state))

        return p
    end
end
```

## Tests

### If production is higher in shallower water 

And reversed with pelagic production.

```{.julia #production-spec}
@testset "Components/Production" begin
    let prod = BenthicProduction(
            maximum_growth_rate = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        input = Input(
            box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
            time = TimeProperties(Δt=1.0u"kyr", steps=10),
            sea_level = t -> 0.0u"m",
            initial_topography = (x, y) -> -10u"m",
            subsidence_rate = 0.0u"m/Myr",
            facies = [Facies(production=prod)],
            insolation = 400.0u"W/m^2")

        state = initial_state(input)
        prod = uniform_production(input)(state)
        @test all(prod[1:end-1,:] .>= prod[2:end,:])
    end

    let prod = PelagicProduction(
            maximum_growth_rate = 5u"1/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        input = Input(
            box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
            time = TimeProperties(Δt=1.0u"kyr", steps=10),
            sea_level = t -> 0.0u"m",
            initial_topography = (x, y) -> -10u"m",
            subsidence_rate = 0.0u"m/Myr",
            facies = [Facies(production=prod)],
            insolation = 400.0u"W/m^2")

        state = initial_state(input)
        prod = uniform_production(input)(state)
        @test all(prod[1:end-1,:] .<= prod[2:end,:])
    end
end
```

### Variable insolation

If insolation increases linearly with time, production at t = 10 should be higher than at t = 1.

```{.julia #production-spec}
@testset "Components/Production/variable_insolation" begin
    let prod = BenthicProduction(
            maximum_growth_rate = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        input = Input(
            box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
            time = TimeProperties(Δt=1.0u"kyr", steps=10),
            sea_level = t -> 0.0u"m",
            initial_topography = (x, y) -> -10u"m",
            subsidence_rate = 0.0u"m/Myr",
            facies = [Facies(production=prod)],
            insolation = t -> 40.0u"W/m^2/kyr" * t)

        state = initial_state(input)
        state.step = 1
        prod1 = copy(uniform_production(input)(state))
        state.step = 10
        prod2 = copy(uniform_production(input)(state))
        @test all(prod2 .> prod1)
    end
end
```

``` {.julia file=test/Components/ProductionSpec.jl}
module ProductionSpec
    using Test
    using CarboKitten
    using CarboKitten.Components.Common
    using CarboKitten.Components.Production: Facies, Input, uniform_production
    using CarboKitten.Components.WaterDepth: initial_state

    <<production-spec>>
end
```
