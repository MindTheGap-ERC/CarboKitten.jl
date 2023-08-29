# ~/~ begin <<docs/src/carbocat.md#src/Burgess2013/Production.jl>>[init]
module Production

using ..Config: Species, Iₖ, k, gₘ

# ~/~ begin <<docs/src/bosscher-1992.md#carbonate-production>>[init]
g(gₘ, I₀, Iₖ, k, w) = gₘ * tanh(I₀/Iₖ * exp(-w * k))
# ~/~ end

function production_rate(I₀::Float64, s::Species, w::Float64)
    w > 0.0 ? g(gₘ(s), I₀, Iₖ(s), k(s), w) : 0.0
end

function production_rate(I₀::Float64, specs::Vector{Species}, spec_map::Matrix{Int}, w::Matrix{Float64})
    production_rate.(I₀, Iterators.map(i -> specs[i], spec_map), w)
end

end
# ~/~ end