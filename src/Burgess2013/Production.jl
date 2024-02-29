# ~/~ begin <<docs/src/carbocat.md#src/Burgess2013/Production.jl>>[init]
module Production

export production_rate

using ..Config: Facies

# ~/~ begin <<docs/src/bosscher-1992.md#carbonate-production>>[init]
g(gₘ, I₀, Iₖ, k, w) = gₘ * tanh(I₀/Iₖ * exp(-w * k))
# ~/~ end

function production_rate(insolation::Float64, facies::F, water_depth::Float64) where {F}
    gₘ = facies.maximum_growth_rate
    I₀ = insolation
    Iₖ = facies.saturation_intensity
    w = water_depth
    k = facies.extinction_coefficient
    return w > 0.0 ? gₘ * tanh(I₀/Iₖ * exp(-w * k)) : 0.0
end

end
# ~/~ end