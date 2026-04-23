# ~/~ begin <<docs/src/active-layer-transport.md#src/Transport/ActiveLayer.jl>>[init]
module ActiveLayer

using Unitful
using StaticArrays
using ...BoundaryTrait
using ...Boxes: Box
using ...Stencil: stencil!

const Rate = typeof(1.0u"m/Myr")
const Amount = typeof(1.0u"m")

function pde_stencil(box::Box{BT}, Δt, ν, out, η, C) where {BT<:Boundary{2}}
    Δx = box.phys_scale
    d = ν * Δt

    stencil!(BT, Size(3, 3), out, η, C) do η, C
        adv = d * ((η[3, 2] - η[1, 2]) * (C[3, 2] - C[1, 2]) +
                   (η[2, 3] - η[2, 1]) * (C[2, 3] - C[2, 1])) /
              (2Δx)^2

        dif = d * C[2, 2] * (η[3, 2] + η[2, 3] + η[1, 2] +
                             η[2, 1] - 4 * η[2, 2]) / (Δx)^2

        prd = C[2, 2]

        max(0.0u"m", adv + dif + prd)
    end
end

end
# ~/~ end