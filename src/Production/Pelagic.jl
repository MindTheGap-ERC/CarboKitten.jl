# ~/~ begin <<docs/src/production/pelagic.md#src/Production/Pelagic.jl>>[init]
module Pelagic

using Unitful
using Interpolations: linear_interpolation
import ..Abstract: AbstractProduction, is_pelagic, production_profile, insolation_curve

@kwdef struct PelagicProduction <: AbstractProduction
    maximum_growth_rate::typeof(1.0u"1/Myr") = 0.0u"1/Myr"
    extinction_coefficient::typeof(1.0u"m^-1") = 0.0u"m^-1"
    saturation_intensity::typeof(1.0u"W/m^2") = 1.0u"W/m^2"
    maximum_production_depth::typeof(1.0u"m") = 200.0u"m"
    table_size::Tuple{Int, Int} = (1000, 1000)
end

is_pelagic(::PelagicProduction) = true

production_profile(input::AbstractInput, p::PelagicProduction) = 
    pelagic_production_lookup(input, p)

end
# ~/~ end
