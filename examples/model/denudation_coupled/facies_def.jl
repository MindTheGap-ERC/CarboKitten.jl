module DenufaciesDef



using CarboKitten.Components
using CarboKitten.Components.Common
using CarboKitten.Models: WithDenudation as WDn
using Unitful

const FACIES = [
    WDn.Facies(
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=50u"m/yr",
        reactive_surface=1000u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.00023u"m/yr"
        ),
    WDn.Facies(
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=25u"m/yr",
        reactive_surface=1000u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.00023u"m/yr"
        ),
    WDn.Facies(
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=12.5u"m/yr",
        reactive_surface=1000u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.00023u"m/yr"
        )
]

export FACIES

end