# ~/~ begin <<docs/src/carbocat.md#src/Burgess2013/Config.jl>>[init]
module Config

using Unitful
using Parameters
export Facies, MODEL1

 struct Facies
    viability_range::Tuple{Int, Int}
    activation_range::Tuple{Int, Int}

    maximum_growth_rate::typeof(1.0u"m/Myr")
    extinction_coefficient::typeof(1.0u"m^-1")
    saturation_intensity::typeof(1.0u"W/m^2")

    reactive_surface::Float64 #reactive surface
    mass_density::Float64 #density of different carb factory
    infiltration_coefficient::Float64 #infiltration coeff
end

const MODEL1 = [
    Facies((4, 10),
           (6, 10),
           500u"m/Myr",
           0.8u"m^-1",
           60u"W/m^2",
           1000,
           2730,
           0.5),

    Facies((4, 10),
           (6, 10),
           400u"m/Myr",
           0.1u"m^-1",
           60u"W/m^2",
           1000,
           2730,
           0.5),

    Facies((4, 10),
           (6, 10),
           100u"m/Myr",
           0.005u"m^-1",
           60u"W/m^2",
           1000,
           2730,
           0.5)
]

end
# ~/~ end