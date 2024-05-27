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

const DenudationIstance = Dissolution(273.0,1000.0,10^(-1.5),2e-3)

const input = Input(
    Box{Shelf}((100, 50), 1.0u"km"),
    TimeProperties(
      1.0u"kyr",
      100,
      1
    ),
    t -> 4.0u"m" * sin(2π * t / 200.0u"kyr"), 
    50.0u"m/Myr",
    p -> x / 300.0,
    MODEL1,
    400.0u"W/m^2",
    DenudationIstance
  )

CaProdErosion.main(input, "data/caps-test.h5")
# ~/~ end