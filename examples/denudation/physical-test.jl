# ~/~ begin <<docs/src/denudation/physical_erosion.md#examples/denudation/physical-test.jl>>[init]
#| requires: examples/denudation/physical-test.jl
#| creates: docs/src/_fig/PhysicalSlope.png
#| collect: figures
module PhysicalSpec

using CairoMakie
using CarboKitten.Denudation.EmpiricalDenudationMod: slope_kernel
using CarboKitten.Denudation.PhysicalErosionMod: physical_erosion, redistribution_kernel
using CarboKitten.Stencil: Boundary, Periodic, offset_value, offset_index, stencil
using CarboKitten.BoundaryTrait
import CarboKitten.Config: Box

const inf = 0.5
const erodibility = 0.0025

function main()
    SLOPE = collect(1:0.2:30)
    DEN_MASS = Vector{Float64}(undef,size(SLOPE))

    for i in eachindex(SLOPE)
        DEN_MASS[i] = physical_erosion(SLOPE[i],inf,erodibility)
    end

    w = 10 .* [ 0.0663001  0.115606  0.646196
    0.601523   0.130196  0.390821
    0.864734   0.902935  0.670354]

    fig = Figure()
    ax = Axis(fig[1,1],xlabel="Slope (degrees)", ylabel="Denudation rates (m/Myr)")
    lines!(ax,SLOPE,DEN_MASS)
    save("docs/src/_fig/PhysicalSlope.png",fig)
end

end

PhysicalSpec.main()
# ~/~ end