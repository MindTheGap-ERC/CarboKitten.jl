# ~/~ begin <<docs/src/production/benthic.md#src/Production/Benthic.jl>>[init]
module Benthic

using Unitful
import .Abstract: is_benthic, insolation_curve, production_profile

# ~/~ begin <<docs/src/production/benthic.md#benthic-production-rate>>[init]
function production_rate(insolation, facies, water_depth)
    gₘ = facies.maximum_growth_rate
    I = insolation / facies.saturation_intensity
    x = water_depth * facies.extinction_coefficient
    return x > 0.0 ? gₘ * tanh(I * exp(-x)) : zero(typeof(gₘ))
end
# ~/~ end
# ~/~ begin <<docs/src/production/benthic.md#benthic-production-rate>>[1]
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

@kwdef struct BenthicProduction <: AbstractProduction
    maximum_growth_rate::typeof(1.0u"m/Myr") = 0.0u"m/Myr"
    extinction_coefficient::typeof(1.0u"m^-1") = 0.0u"m^-1"
    saturation_intensity::typeof(1.0u"W/m^2") = 1.0u"W/m^2"
end

is_benthic(::BenthicProduction) = true

function production_profile(input::AbstractInput, p::BenthicProduction)
    I_of_t = insolation_curve(input)
    return (t, w) -> benthic_production(I_of_t(t), p, w)
end

end
# ~/~ end
