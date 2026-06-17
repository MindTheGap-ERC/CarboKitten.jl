# ~/~ begin <<docs/src/components/production.md#src/Production.jl>>[init]
module Production

using Unitful
using QuadGK
using Interpolations

import ..CarboKitten: AbstractInput, time_axis
using ..Utility: in_units_of

export production_profile

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
and cap by available accommodation. Returns the deposited thickness for `dt`.
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
# ~/~ begin <<docs/src/components/production.md#insolation-curve>>[init]
"""
    insolation_curve(input) -> function(time) -> insolation

Build a closure mapping time to insolation from the input specification.
Handles constant (`Quantity`), tabular (`AbstractVector`), and functional
insolation inputs.
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
# ~/~ begin <<docs/src/components/production.md#production-profile>>[init]
"""
    production_profile(input::AbstractInput, p)

Given an input and a production configuration, returns a function
`(time, water_depth) -> production_rate`.

Insolation is read from `input` and composed into the returned closure via
`insolation_curve` — the caller does not need to pass insolation at each
time step.

The default implementation assumes `p` is already a callable
`(time, water_depth) -> rate`.

## Example

    struct MyProduction <: AbstractProduction
        ...
    end

    import CarboKitten.Production: production_profile

    production_profile(input::AbstractInput, p::MyProduction) =
        function (time, water_depth)
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

"""
    is_interpolated(obj)

Predicate to determine if a facies or production spec is interpolation-based.
Defaults to `false`.
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
production_profile(input::AbstractInput, p::PelagicProduction) = pelagic_production_lookup(input, p)
# ~/~ end
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

`depth_knots` need not be sorted; they are sorted internally.

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

is_benthic(::InterpolatedProduction)      = false
is_pelagic(::InterpolatedProduction)      = false
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

Wraps `base::AbstractProduction`, multiplying its output by `factor` during
`t_range`. Outside `t_range` the base production is unchanged.

This implements the modifier pattern as `AbstractProduction -> AbstractProduction`:
modifiers compose directly in the production spec rather than in a separate
`production_modifiers` list on `Input`.
"""
@kwdef struct MultiplyProduction <: AbstractProduction
    base::AbstractProduction
    factor::Float64
    t_range::_ProdTimeSpec = (:)
end
MultiplyProduction(base, factor::Real; kwargs...) =
    MultiplyProduction(; base=base, factor=Float64(factor), kwargs...)

is_benthic(p::MultiplyProduction)      = is_benthic(p.base)
is_pelagic(p::MultiplyProduction)      = is_pelagic(p.base)
is_interpolated(p::MultiplyProduction) = is_interpolated(p.base)

function production_profile(input::AbstractInput, p::MultiplyProduction)
    base_profile = production_profile(input, p.base)
    return function(t, w)
        f = p.t_range isa Colon || (p.t_range[1] <= t <= p.t_range[2]) ? p.factor : 1.0
        return base_profile(t, w) * f
    end
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
