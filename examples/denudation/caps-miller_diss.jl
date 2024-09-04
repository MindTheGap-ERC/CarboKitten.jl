# ~/~ begin <<docs/src/ca-with-production.md#examples/caps-osc.jl>>[init]
using CarboKitten
using CarboKitten.Model.WithDenudation
using CarboKitten.Burgess2013
using CarboKitten.Model.WithDenudation: Input
using CSV
using DataFrames
using Interpolations
using Unitful
using CarboKitten.Config: Box, Vectors, TimeProperties
using CarboKitten.BoundaryTrait
using CarboKitten.Denudation: DissolutionMod, NoDenudationMod, PhysicalErosionMod, EmpiricalDenudationMod
using CarboKitten.Denudation: Dissolution, EmpiricalDenudation, PhysicalErosion, NoDenudation

const MODEL1 = [
    WithDenudation.Facies(viability_range = (4, 10),
    activation_range = (6, 10),
    maximum_growth_rate = 500u"m/Myr",
    extinction_coefficient = 0.8u"m^-1",
    saturation_intensity = 60u"W/m^2",
    reactive_surface = 1000u"m^2/m^3",
    mass_density = 2730u"kg/m^3",
    infiltration_coefficient= 0.5),

    WithDenudation.Facies(viability_range = (4, 10),
    activation_range = (6, 10),
    maximum_growth_rate = 400u"m/Myr",
    extinction_coefficient = 0.1u"m^-1",
    saturation_intensity = 60u"W/m^2",
    reactive_surface = 1000u"m^2/m^3",
    mass_density = 2730u"kg/m^3",
    infiltration_coefficient= 0.5),

    WithDenudation.Facies(viability_range = (4, 10),
    activation_range = (6, 10),
    maximum_growth_rate = 100u"m/Myr",
    extinction_coefficient = 0.005u"m^-1",
    saturation_intensity = 60u"W/m^2",
    reactive_surface = 1000u"m^2/m^3",
    mass_density = 2730u"kg/m^3",
    infiltration_coefficient= 0.5)
]

const DENUDATION = Dissolution(temp = 273.0u"K",precip = 1.0u"m", pco2 = 10^(-1.5)*u"atm",reactionrate = 2e-3u"m/yr")#PhysicalErosionParam(1e-3)EmpericalDenudationParam(1.0)

const INPUT = Input(
    box=Box{Shelf}(grid_size = (100, 50), phys_scale = 1.0u"km"),
    time=TimeProperties(
      Δt = 1.0u"kyr",
      steps = 100,
      write_interval = 1
    ),
    sea_level = t -> 0.0, 
    subsidence_rate = 50.0u"m/Myr",
    initial_depth = x -> x/300.0,
    facies = MODEL1,
    insolation = 400.0u"W/m^2",
    denudation = DENUDATION
  )

WithDenudation.main(INPUT, "data/caps-test2.h5")
# ~/~ end