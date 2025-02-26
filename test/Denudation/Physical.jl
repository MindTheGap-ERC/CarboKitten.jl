module PhysicalUnitTest
using Test
using CarboKitten.Denudation.PhysicalErosionMod: redistribution_kernel, mass_erosion, total_mass_redistribution
using CarboKitten: Box
using CarboKitten.Stencil: Periodic, Reflected, stencil
using Unitful

WD = [ 0.0663001  0.115606  0.646196
0.601523   0.130196  0.390821
0.864734   0.902935  0.670354]

WD_zeros = zeros(3,3)

Redis_kernel_result = Array{Float64,2}(undef,3,3)
Redis_kernel_result_zeros = Array{Float64,2}(undef,3,3)

Redis_kernel_result = redistribution_kernel(WD,1.0)
Redis_kernel_result_zeros = redistribution_kernel(WD_zeros,1.0)

# test the normalization of the redist. kernel
@test sum(Redis_kernel_result) == 1.0
# test all zero waterdepth
@test Redis_kernel_result_zeros == zeros(3,3)


const BOX = Box{Periodic{2}}(grid_size=(3, 3), phys_scale=1.0u"km")
const DENU_MASS = 2*ones(4,4)

# the mass redistributed around the cell should be the same as the mass loss in the cell
@test mass_erosion(BOX, DENU_MASS, WD, CartesianIndex(2,2)) |> sum ≈ DENU_MASS[2,2] 

# test the conservation of mass during redistribution, on a bigger grid
const BOX_T = Box{Periodic{2}}(grid_size=(4, 4), phys_scale=1.0u"km")
MASS = 100*rand(Float64, 4,4)#Array{Float64,2}(undef,3,3)
MASS_PRE = deepcopy(MASS)
WD_T = [ 0.0663001  0.115606  0.646196 0.45
0.601523   0.130196  0.390821 0.56
0.864734   0.902935  0.670354 0.393
0.564734   0.102935  0.473354 0.893]

# the total mass on grid + the pre-removed denudated mass before 
# should match the total mass after, once the denudated mass is redistributed
# this currently fails - local depressions in the slope field break conservation of mass!
total_mass_redistribution(BOX_T, DENU_MASS, WD_T, MASS)
@test sum(MASS) ≈ sum(MASS_PRE) + sum(DENU_MASS) # - DENU_MASS[2,2]

end