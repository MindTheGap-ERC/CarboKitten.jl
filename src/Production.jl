# ~/~ begin <<docs/src/production/production.md#src/Production.jl>>[init]
module Production

include("Production/Abstract.jl")
include("Production/Insolation.jl")
include("Production/Benthic.jl")
include("Production/Pelagic.jl")
include("Production/Interpolated.jl")
include("Production/Modifiers.jl")

using Unitful

import .Abstract: AbstractProduction, NoProduction, production_profile, is_benthic, is_pelagic, is_interpolated
import .Benthic: BenthicProduction
import .Pelagic: PelagicProduction
import .Interpolated: InterpolatedProduction
import .Modifiers: MultiplyProduction, ProductionBoost
import .Insolation: insolation_curve

export AbstractProduction, production_profile, NoProduction, BenthicProduction, PelagicProduction,
    InterpolatedProduction, MultiplyProduction, ProductionBoost, EXAMPLE

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
        saturation_intensity=60u"W/m^2"),
    :interpolated => InterpolatedProduction(
        maximum_production=500u"m/Myr",
        depth_knots=[0.0u"m", 5.0u"m", 15.0u"m", 30.0u"m", 50.0u"m"],
        multipliers=[0.0, 1.0, 1.0, 0.4, 0.0]),
    :time_varying => MultiplyProduction(
        BenthicProduction(
            maximum_growth_rate=500u"m/Myr",
            extinction_coefficient=0.8u"m^-1",
            saturation_intensity=60u"W/m^2"),
        0.5;
        t_range=(0.5u"Myr", 1.0u"Myr"))
)

end
# ~/~ end
