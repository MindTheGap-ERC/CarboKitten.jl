# ~/~ begin <<docs/src/carbocat.md#src/Burgess2013/Config.jl>>[init]
module Config

using Unitful
export Facies, MODEL1

@kwdef struct Facies
    viability_range::Tuple{Int, Int}
    activation_range::Tuple{Int, Int}

    maximum_growth_rate::typeof(1.0u"m/Myr")
    extinction_coefficient::typeof(1.0u"m^-1")
    saturation_intensity::typeof(1.0u"W/m^2")

end

# const MODEL1 = [
#     Facies(viability_range = (4, 10),
#            activation_range = (6, 10),
#            maximum_growth_rate = 500u"m/Myr",
#            extinction_coefficient = 0.8u"m^-1",
#            saturation_intensity = 60u"W/m^2"),

#     Facies(viability_range = (4, 10),
#            activation_range = (6, 10),
#            maximum_growth_rate = 400u"m/Myr",
#            extinction_coefficient = 0.1u"m^-1",
#            saturation_intensity = 60u"W/m^2"),

#     Facies(viability_range = (4, 10),
#            activation_range = (6, 10),
#            maximum_growth_rate = 100u"m/Myr",
#            extinction_coefficient = 0.005u"m^-1",
#            saturation_intensity = 60u"W/m^2")
# ]

end
# ~/~ end
