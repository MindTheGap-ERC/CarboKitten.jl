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

function capped_production(f, insolation, water_depth, dt)
    clip_positive(x::T) where {T} = max(x, zero(T))
    p = clip_positive(f(insolation, water_depth))
    return min(max(0.0u"m", water_depth), p * dt)
end
# ~/~ end
# ~/~ begin <<docs/src/components/production.md#pelagic-production>>[init]
function pelagic_production(insolation, facies, water_depth)
    return quadgk(w -> production_rate(insolation, facies, w), 0.0u"m", water_depth)[1]
end
# ~/~ end
# ~/~ begin <<docs/src/components/production.md#production-profile>>[init]
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
function is_benthic end

"""
    is_pelagic(obj)

Predicate to determine if a facies or production spec is pelagic.
Defaults to `false`.
"""
function is_pelagic end

abstract type AbstractProduction end

struct NoProduction <: AbstractProduction
end

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

# ~/~ end
# ~/~ begin <<docs/src/components/production.md#production-lookup>>[init]
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
# ~/~ end

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
        saturation_intensity=60u"W/m^2")
)

end
# ~/~ end
