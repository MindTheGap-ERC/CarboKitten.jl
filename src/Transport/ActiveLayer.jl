# ~/~ begin <<docs/src/active-layer-transport.md#src/Transport/ActiveLayer.jl>>[init]
module ActiveLayer

using Unitful
using ...BoundaryTrait
using ...Boxes: Box
using ...Stencil: stencil

const Rate = typeof(1.0u"m/Myr")
const Amount = typeof(1.0u"m")

function pde_stencil(box::Box{BT}, Δt::Unitful.Time{Float64}, ν::Unitful.Velocity{Float64}) where {BT<:Boundary{2}}
    Δx = box.phys_scale
    d = ν * Δt

    function kernel(x)
        adv = d * ((x[3, 2][1] - x[1, 2][1]) * (x[3, 2][2] - x[1, 2][2]) +
                   (x[2, 3][1] - x[2, 1][1]) * (x[2, 3][2] - x[2, 1][2])) /
              (2Δx)^2

        dif = d * x[2, 2][2] * (x[3, 2][1] + x[2, 3][1] + x[1, 2][1] +
                                x[2, 1][1] - 4 * x[2, 2][1]) / (Δx)^2

        prd = x[2, 2][2]

        return max(0.0u"m", adv + dif + prd)
    end

    stencil(Tuple{Amount,Amount}, Amount, BT, (3, 3), kernel)
end

end
# ~/~ end
