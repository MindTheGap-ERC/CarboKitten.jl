# ~/~ begin <<docs/src/components/production.md#src/Production.jl>>[init]
module Production

using Unitful
using QuadGK
using Interpolations

import ..CarboKitten: AbstractInput, time_axis
using ..Utility: in_units_of

export production_profile

# =============================================================================
# Core production equation (Bosscher & Schlager 1992)
# =============================================================================

# ~/~ begin <<docs/src/components/production.md#component-production-rate>>[init]
function production_rate(insolation, facies, water_depth)
    gₘ = facies.maximum_growth_rate
    I = insolation / facies.saturation_intensity
    x = water_depth * facies.extinction_coefficient
    return x > 0.0 ? gₘ * tanh(I * exp(-x)) : zero(typeof(gₘ))
end

benthic_production(i, f, w) = production_rate(i, f, w)

"""
    capped_production(f, time, water_depth, dt)

Apply production function `f(time, water_depth) -> rate`, clip to non-negative,
and cap by available accommodation (water depth). Returns the deposited thickness
for time step `dt`.
"""
function capped_production(f, time, water_depth, dt)
    clip_positive(x::T) where {T} = max(x, zero(T))
    p = clip_positive(f(time, water_depth))
    return min(max(0.0u"m", water_depth), p * dt)
end
# ~/~ end

# ~/~ begin <<docs/src/components/production.md#pelagic-production>>[init]
function pelagic_production(insolation, facies, water_depth)
    return quadgk(w -> production_rate(insolation, facies, w), 0.0u"m", water_depth)[1]
end
# ~/~ end

# =============================================================================
# Insolation helper
# =============================================================================

# ~/~ begin <<docs/src/components/production.md#insolation-curve>>[init]
"""
    insolation_curve(input) -> function(time) -> insolation

Build a closure that maps time to insolation from the input specification.
This is used by `production_profile` to capture insolation inside the
returned `(time, water_depth) -> rate` closure, removing the need to pass
insolation explicitly through the model loop.
"""
function insolation_curve(input::AbstractInput)
    insolation_param = input.insolation
    if insolation_param isa Quantity
        return _ -> insolation_param
    elseif insolation_param isa AbstractVector
        t_axis = time_axis(input)
        t_vals = ustrip.(u"Myr", t_axis)
        I_vals = ustrip.(u"W/m^2", insolation_param)
        itp = linear_interpolation(t_vals, I_vals, extrapolation_bc=Flat())
        return t -> itp(ustrip(u"Myr", t)) * u"W/m^2"
    else
        return t -> insolation_param(t)
    end
end
# ~/~ end

# =============================================================================
# Production profile API
# =============================================================================

# ~/~ begin <<docs/src/components/production.md#production-profile>>[init]
"""
    production_profile(input::AbstractInput, p)

Given an input and a production configuration, returns a function
`(time, water_depth) -> production_rate`.

Insolation is read from `input` and composed into the returned closure —
the caller does not need to pass insolation at each time step.

The default implementation assumes `p` is already a callable
`(time, water_depth) -> rate`.
"""
production_profile(input::AbstractInput, p) = p

"""
    is_benthic(obj)

Predicate: is this a benthic (seafloor) production spec? Default `false`.
"""
is_benthic(p) = false

"""
    is_pelagic(obj)

Predicate: is this a pelagic (water column) production spec? Default `false`.
"""
is_pelagic(p) = false

"""
    is_interpolated(obj)

Predicate: is this a knot-based interpolated production spec? Default `false`.
"""
is_interpolated(p) = false

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

"""
    production_profile(input, p::BenthicProduction)

Returns `(time, water_depth) -> rate`. Insolation at `time` is looked up
from `input.insolation` via `insolation_curve`.
"""
function production_profile(input::AbstractInput, p::BenthicProduction)
    I_of_t = insolation_curve(input)
    return (t, w) -> benthic_production(I_of_t(t), p, w)
end

@kwdef struct PelagicProduction <: AbstractProduction
    maximum_growth_rate::typeof(1.0u"1/Myr") = 0.0u"1/Myr"
    extinction_coefficient::typeof(1.0u"m^-1") = 0.0u"m^-1"
    saturation_intensity::typeof(1.0u"W/m^2") = 1.0u"W/m^2"
    maximum_production_depth::typeof(1.0u"m") = 200.0u"m"
    table_size::Tuple{Int, Int} = (1000, 1000)
end

is_benthic(::PelagicProduction) = false
is_pelagic(::PelagicProduction) = true

function production_profile(input::AbstractInput, p::PelagicProduction)
    return pelagic_production_lookup(input, p)
end
# ~/~ end

# =============================================================================
# Input parameters
# =============================================================================

# ~/~ begin <<docs/src/components/production.md#production-input>>[init]
@kwdef struct Input <: AbstractInput
    insolation
end

@kwdef struct Facies <: AbstractFacies
    production = NoProduction()
end

is_benthic(facies::AbstractFacies) = is_benthic(facies.production)
is_pelagic(facies::AbstractFacies) = is_pelagic(facies.production)
# ~/~ end

# =============================================================================
# Pelagic lookup table
# =============================================================================

# ~/~ begin <<docs/src/components/production.md#production-lookup>>[init]
function pelagic_production_lookup(input::AbstractInput, prod::PelagicProduction)
    I_of_t = insolation_curve(input)
    depth_grid = LinRange(0.0, prod.maximum_production_depth |> in_units_of(u"m"), prod.table_size[2])

    if input.insolation isa Quantity
        # Constant insolation — 1D depth lookup, time argument ignored
        I0 = input.insolation
        production_values = [pelagic_production(I0, prod, w * u"m") |> in_units_of(u"m/Myr")
                             for w in depth_grid]
        itp = linear_interpolation(depth_grid, production_values, extrapolation_bc=Line())
        return (_, w) -> itp(w |> in_units_of(u"m")) * u"m/Myr"
    end

    # Variable insolation — 2D (insolation × depth) lookup
    t_axis = time_axis(input)
    I_vals = [I_of_t(t) |> in_units_of(u"W/m^2") for t in t_axis]
    I_min, I_max = extrema(I_vals)

    insolation_grid = LinRange(I_min, I_max, prod.table_size[1])
    production_values = [
        pelagic_production(I * u"W/m^2", prod, w * u"m") |> in_units_of(u"m/Myr")
        for I in insolation_grid, w in depth_grid
    ]
    itp = linear_interpolation((collect(insolation_grid), collect(depth_grid)),
                               production_values, extrapolation_bc=Line())
    return (t, w) -> itp(I_of_t(t) |> in_units_of(u"W/m^2"), w |> in_units_of(u"m")) * u"m/Myr"
end
# ~/~ end

# =============================================================================
# Interpolated (knot-based) production curve
# =============================================================================

"""
    InterpolatedProduction(; maximum_production, depth_knots, multipliers)

A depth-only production curve defined by a peak rate and a piecewise-linear
shape over `(depth, multiplier)` knots. Independent of insolation.

    rate(t, w) = maximum_production × interpolate(depth_knots, multipliers; w)

# Example

    InterpolatedProduction(
        maximum_production = 500.0u"m/Myr",
        depth_knots        = [0.0u"m", 5.0u"m", 20.0u"m", 50.0u"m"],
        multipliers        = [0.0,     1.0,     0.6,      0.0])
"""
@kwdef struct InterpolatedProduction <: AbstractProduction
    maximum_production::typeof(1.0u"m/Myr") = 0.0u"m/Myr"
    depth_knots::Vector{typeof(1.0u"m")}    = typeof(1.0u"m")[]
    multipliers::Vector{Float64}            = Float64[]
end

is_benthic(::InterpolatedProduction)     = false
is_pelagic(::InterpolatedProduction)     = false
is_interpolated(::InterpolatedProduction) = true

function production_profile(::AbstractInput, p::InterpolatedProduction)
    @assert length(p.depth_knots) == length(p.multipliers)
    @assert length(p.depth_knots) >= 2
    depths_m = [d |> in_units_of(u"m") for d in p.depth_knots]
    order = sortperm(depths_m)
    itp = linear_interpolation(depths_m[order], p.multipliers[order], extrapolation_bc=Flat())
    max_rate = p.maximum_production
    return (_, w) -> max_rate * itp(w |> in_units_of(u"m"))
end

# =============================================================================
# Time-window modifier — AbstractProduction transformer
# =============================================================================

const _ProdTime     = typeof(1.0u"Myr")
const _ProdTimeSpec = Union{Colon, Tuple{_ProdTime,_ProdTime}}

"""
    MultiplyProduction(base, factor; t_range=:)

Transforms `base::AbstractProduction` into a new production spec that
multiplies the output by `factor` during `t_range`. Outside `t_range` the
base production is returned unchanged.

This implements the modifier pattern as an `AbstractProduction -> AbstractProduction`
transformation: modifiers compose directly in the production spec rather than
living in a separate `production_modifiers` list.

# Example

    # Reef growth halved between 0.5 and 1.0 Myr
    MultiplyProduction(
        BenthicProduction(maximum_growth_rate=500u"m/Myr", ...),
        0.5;
        t_range=(0.5u"Myr", 1.0u"Myr"))

    # Composing multiple modifiers
    MultiplyProduction(
        MultiplyProduction(base_prod, 0.5; t_range=(0u"Myr", 0.3u"Myr")),
        2.0;
        t_range=(0.8u"Myr", 1.0u"Myr"))
"""
@kwdef struct MultiplyProduction <: AbstractProduction
    base::AbstractProduction
    factor::Float64
    t_range::_ProdTimeSpec = (:)
end
MultiplyProduction(base, factor::Real; kwargs...) =
    MultiplyProduction(; base=base, factor=Float64(factor), kwargs...)

# Delegate predicates to the wrapped base production
is_benthic(p::MultiplyProduction)     = is_benthic(p.base)
is_pelagic(p::MultiplyProduction)     = is_pelagic(p.base)
is_interpolated(p::MultiplyProduction) = is_interpolated(p.base)

function production_profile(input::AbstractInput, p::MultiplyProduction)
    base_profile = production_profile(input, p.base)
    return function(t, w)
        f = p.t_range isa Colon || (p.t_range[1] <= t <= p.t_range[2]) ? p.factor : 1.0
        return base_profile(t, w) * f
    end
end

# =============================================================================
# Example configurations
# =============================================================================

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
        multipliers=[0.0, 1.0, 1.0, 0.4, 0.0]),
    :time_varying => MultiplyProduction(
        BenthicProduction(
            maximum_growth_rate=500u"m/Myr",
            extinction_coefficient=0.8u"m^-1",
            saturation_intensity=60u"W/m^2"),
        0.5;
        t_range=(0.5u"Myr", 1.0u"Myr"))
)

end
# ~/~ end
