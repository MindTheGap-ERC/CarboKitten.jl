# ~/~ begin <<docs/src/carbocat.md#src/Burgess2013/Production.jl>>[init]
module Production

export production_rate

using Unitful
using ..Config: Facies

# ~/~ begin <<docs/src/bosscher-1992.md#carbonate-production>>[init]
g(gₘ, I₀, Iₖ, k, w) = gₘ * tanh(I₀/Iₖ * exp(-w * k))
# ~/~ end

function production_rate(insolation, facies, water_depth)
    gₘ = facies.maximum_growth_rate
    I = insolation / facies.saturation_intensity
    x = water_depth * facies.extinction_coefficient
    return water_depth > 0.0u"m" ? gₘ * tanh(I * exp(-x)) : 0.0u"m/Myr"
end

end
# ~/~ end