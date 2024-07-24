module WithDenudation
using Unitful

@kwdef struct Facies
    viability_range::Tuple{Int,Int}
    activation_range::Tuple{Int,Int}

    maximum_growth_rate::typeof(1.0u"m/Myr")
    extinction_coefficient::typeof(1.0u"m^-1")
    saturation_intensity::typeof(1.0u"W/m^2")

    reactive_surface::Float64 #reactive surface
    mass_density::Float64 #density of different carb factory
    infiltration_coefficient::Float64 #infiltration coeff
end

end
