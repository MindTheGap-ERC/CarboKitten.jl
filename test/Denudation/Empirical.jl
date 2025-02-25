module EmpiricalUnitTest
using Test
using CarboKitten.Denudation.EmpiricalDenudationMod: slope_kernel

w = 10 .* ones(3,3)
w1 = rand(3,3)
w2 = 10 * rand(3,3)
w3 = ones(3,3)
w3[5] = 0
@test slope_kernel(w,1.0) == 0.0
@test slope_kernel(w1,1.0) < slope_kernel(w2,1.0) 
@test slope_kernel(w3,1.0) == 0.0

end