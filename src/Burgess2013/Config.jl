# ~/~ begin <<docs/src/carbocat.md#src/Burgess2013/Config.jl>>[init]
module Config

struct Species
    viability_range::Tuple{Int, Int}
    activation_range::Tuple{Int, Int}

    maximum_growth_rate::Float64
    extinction_coefficient::Float64
    saturation_intensity::Float64
end

Iₖ(s::Species) = s.saturation_intensity
k(s::Species) = s.extinction_coefficient
gₘ(s::Species) = s.maximum_growth_rate

model1 = [
    Species((4, 10), (6, 10), 500, 0.8, 300),
    Species((4, 10), (6, 10), 400, 0.1, 300),
    Species((4, 10), (6, 10), 100, 0.005, 300)
]

end
# ~/~ end