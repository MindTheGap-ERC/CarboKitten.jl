module PhysicalUnitTest
using Test
using CarboKitten.Denudation.PhysicalErosionMod: redistribution_kernel, mass_erosion, total_mass_redistribution
using CarboKitten: Box
using CarboKitten.Stencil: Periodic, Reflected, stencil
using CarboKitten.Denudation.EmpiricalDenudationMod: slope_kernel

WD = [ 0.0663001  0.115606  0.646196
0.601523   0.130196  0.390821
0.864734   0.902935  0.670354]

WD_zeros = zeros(3,3)

Redis_kernel_result = Array{Float64,2}(undef,3,3)
Redis_kernel_result_zeros = Array{Float64,2}(undef,3,3)

Redis_kernel_result = redistribution_kernel(WD,1.0)
Redis_kernel_result_zeros = redistribution_kernel(WD_zeros,1.0)

SLOPE = rand(Float64, 3,3)
MASS = Array{Float64,2}(undef,3,3)

@test sum(Redis_kernel_result) == 1.0
@test Redis_kernel_result_zeros == zeros(3,3)

const BOX = Box{Periodic{2}}(grid_size=(3, 3), phys_scale=1.0u"km")
const DENU_MASS = ones(3,3)

slopefn = stencil(Float64, Periodic{2}, (3, 3), slope_kernel)
slopefn(WD, SLOPE, BOX.phys_scale ./u"m")

@test sum(SLOPE) > 0.0

@test mass_erosion(BOX, DENU_MASS, WD, CartesianIndex(2,2)) |> sum > 0.0
@test mass_erosion(BOX, DENU_MASS, WD, CartesianIndex(2,2)) |> sum ≈ DENU_MASS[2,2] 

total_mass_redistribution(BOX, DENU_MASS, WD, MASS)
@test MASS |> sum ≈ sum(DENU_MASS) - DENU_MASS[2,2]

end