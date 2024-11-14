# ~/~ begin <<docs/src/onshore-transport.md#src/Components/ActiveLayerOnshore.jl>>[init]
@compose module ActiveLayerOnshore
@mixin WaterDepth, FaciesBase, SedimentBuffer
using ..Common

struct Facies <: AbstractFacies
    flux_differential
end

function pde_stencil(box::Box{BT}, ν, w) where {BT <: Boundary{2}}
    Δx = box.phys_scale

    function kernel(x)
        adv = ν * ((x[3, 2][1] - x[1, 2][1]) * ((x[3, 2][2] - x[1, 2][2]) - w[1]) +
                   (x[2, 3][1] - x[2, 1][1]) * ((x[2, 3][2] - x[2, 1][2]) - w[2])) /
                  (2Δx)^2

        dif = ν * x[2, 2][2] * (x[3, 2][1] + x[2, 3][1] + x[1, 2][1] +
                  x[2, 1][1] - 4*x[2, 2][1]) / (Δx)^2

        prd = x[2, 2][2]

        return max(0.0u"m", adv + dif + prd)
    end

    stencil(Tuple{Amount, Amount}, Amount, BT, (3, 3), kernel)
end
end
# ~/~ end
