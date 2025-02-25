module DissolutionUnitTest
using Test
using CarboKitten: Box
using CarboKitten.Stencil: Periodic, Reflected, stencil
using CarboKitten.Denudation.EmpiricalDenudationMod: slope_kernel
using CarboKitten.Denudation.DissolutionMod: denudation, Dissolution
using Unitful

WD = [0.0663001  0.115606  0.646196
0.601523   0.130196  0.390821
0.864734   0.902935  0.670354]

@kwdef struct facies
    infiltration_coefficient :: Float64
    mass_density :: typeof(1.0u"kg/m^3")
    reactive_surface :: typeof(1.0u"m^2/m^3")
end

@kwdef struct state
    ca::Matrix{Int}
end

const DIS = Dissolution(temp = 285u"K",
precip  = 1.0u"m",
pco2 = 10^(-1.5)*1.0u"atm",
reactionrate = 0.1u"m/yr")

const STATE = state(ca = 
[ 0  0  0
0  1  0 
0  0  1 ])

const Facies = [facies(infiltration_coefficient = 0.5,
mass_density = 2.73u"kg/m^3",
reactive_surface = 1000u"m^2/m^3")]

SLOPE = rand(Float64, 3,3)
BOX = Box{Periodic{2}}(grid_size=(3, 3), phys_scale=1.0u"km")
slopefn = stencil(Float64, Periodic{2}, (3, 3), slope_kernel)
slopefn(WD, SLOPE, BOX.phys_scale ./u"m")
denudation(BOX, DIS, WD, SLOPE, Facies, STATE)
@test (denudation(BOX, DIS, WD, SLOPE, Facies, STATE) |> sum) ≈ 0.0u"m/yr"

end