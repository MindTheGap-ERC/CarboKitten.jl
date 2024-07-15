# ~/~ begin <<docs/src/ca-with-production.md#examples/caps-osc.jl>>[init]
using CarboKitten
using CarboKitten.CaProdErosion
using CarboKitten.Burgess2013
using CarboKitten.InputConfig: Input
using CSV
using DataFrames
using Interpolations
using Unitful
using CarboKitten.Config: Box, Vectors, TimeProperties
using CarboKitten.BoundaryTrait
using CarboKitten.Denudation: Dissolution, NoDenudation, PhysicalErosionParam, EmpericalDenudationParam

function sealevel_curve(t,filepath)
    data = DataFrame(CSV.File(filepath))
    data = hcat(collect(0:length(data[:,1]).-1)./1000, data[:,1]) #The output sealevel curve from R does not have time tab and have to add it 
    x = linear_interpolation(data[:,1], data[:,2])
    return x(t)
end 

const MODEL1 = [
    Burgess2013.Facies(viability_range = (4, 10),
    activation_range = (6, 10),
    maximum_growth_rate = 500u"m/Myr",
    extinction_coefficient = 0.8u"m^-1",
    saturation_intensity = 60u"W/m^2",
    reactive_surface = 1000,
    mass_density = 2730,
    infiltration_coefficient= 0.5),

    Burgess2013.Facies(viability_range = (4, 10),
    activation_range = (6, 10),
    maximum_growth_rate = 400u"m/Myr",
    extinction_coefficient = 0.1u"m^-1",
    saturation_intensity = 60u"W/m^2",
    reactive_surface = 1000,
    mass_density = 2730,
    infiltration_coefficient= 0.5),

    Burgess2013.Facies(viability_range = (4, 10),
    activation_range = (6, 10),
    maximum_growth_rate = 100u"m/Myr",
    extinction_coefficient = 0.005u"m^-1",
    saturation_intensity = 60u"W/m^2",
    reactive_surface = 1000,
    mass_density = 2730,
    infiltration_coefficient= 0.5)
]

const DENUDATION = EmpericalDenudationParam(1.0)#PhysicalErosionParam(1e-3)

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
    denudationparam = DENUDATION
  )

CaProdErosion.main(INPUT, "data/caps-test.h5")
# ~/~ end