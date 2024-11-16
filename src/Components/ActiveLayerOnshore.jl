# ~/~ begin <<docs/src/onshore-transport.md#src/Components/ActiveLayerOnshore.jl>>[init]
@compose module ActiveLayerOnshore
@mixin WaterDepth, FaciesBase, SedimentBuffer, ActiveLayer
using ..Common
using ...Stencil: stencil

# ~/~ begin <<docs/src/onshore-transport.md#onshore-facies>>[init]
struct Facies <: AbstractFacies
    onshore_velocity
end
# ~/~ end

"""
    pde_stencil(box, ν, wf)

Creates a stencil function for active layer onshore transport, where
`d` is the diffusion coefficient and `wf` is a function of waterdepth
returning the sediment velocity and its derivative (shear) as 2-vectors.

The resulting stencil acts on an array of `Tuple{Length, Amount}`, being
waterdepth and entrained sediment, writing to an array of `Amount` being
deposited sediment.
"""
function onshore_transport_stencil(box::Box{BT}, ν, Δt, sf::F) where {BT<:Boundary{2},F}
    Δx = box.phys_scale

    function kernel(x)
        (w, P) = x[2, 2][2]
        sv, ss = sf(w)
        d = ν * Δt

        adv = -((x[3, 2][1] - x[1, 2][1]) / (2Δx) * (d * (x[3, 2][2] - x[1, 2][2]) / (2Δx) - P * ss[1]) +
                (x[2, 3][1] - x[2, 1][1]) / (2Δx) * (d * (x[2, 3][2] - x[2, 1][2]) / (2Δx) - P * ss[2]))

        dif = -d * P * (x[3, 2][1] + x[2, 3][1] + x[1, 2][1] +
                        x[2, 1][1] - 4 * x[2, 2][1]) / (Δx)^2

        prd = sv[1] * (x[3, 2] - x[1, 2]) / (2Δx) + sv[2] * (x[2, 3] - x[2, 1]) / (2Δx) + x[2, 2][2]

        return max(0.0u"m", adv + dif + prd)
    end

    stencil(Tuple{Amount,Amount}, Amount, BT, (3, 3), kernel)
end

"""
    transportation(input)

Computes the transport using Active Layer with Onshore vector.
"""
function transportation(input)
    w = water_depth(input)

    # We always return this array
    transported_output = Array{Amount,3}(undef, n_facies(input), input.box.grid_size...)
    stencils = [
        let stc = onshore_transport_stencil(input.box, f.diffusion_coefficient, f.onshore_velocity)
            (w, p) -> @views stc(tuple.(w, p[i, :, :]), transported_output[i, :, :])
        end for (i, f) in enumerate(input.facies)]

    return function (state, active_layer::Array{Amount,3})
        wd = w(state)

        for stc in stencils
            stc(wd, active_layer)
        end

        return transported_output
    end
end

end
# ~/~ end
