# Onshore transport

In the Active Layer approach, we derived a diffusion equation for our transport model starting with a proposition for the sediment flux as,

$$q_f = -\nu_f P_f \nabla \eta,$$

where $\nu_f$ is the diffusion coefficient, $P_f$ the production and $\eta$ the sediment height.
Then inserting that into the mass balance equation,

$$\partial_t \eta_f = -\nabla \cdot q_f + P_f.$$

Now we'll add a vector component $P_f w_f(\eta)$ to the flux:

$$q_f = -\nu_f P_f \nabla \eta + P_f w_f(\eta).$$

This leads to the modified advection-diffusion equation:

$$\partial_t \eta_f = \nu_f (P_f \nabla^2 \eta + \nabla P_f \cdot \nabla \eta) - \nabla \cdot (P_f w_f(\eta)) + P_f.$$

Using the chain-rule $\nabla \cdot (P_f w_f(\eta)) = P_f (w_f' \cdot \nabla \eta) + \nabla P_f \cdot w_f$, we arrive at:

$$\partial_t \eta_f = \nu_f P_f \nabla^2 \eta + (\nu_f \nabla P_f - P_f w_f') \cdot \nabla \eta - \nabla P_f \cdot w_f + P_f.$$

So we modify the advection component in the active layer approach with $P_f w_f'$, the derivative of the wave induced flux with respect to water depth and add a new term $\nabla P_f \cdot w_f$.

``` {.julia file=src/Components/ActiveLayerOnshore.jl}
@compose module ActiveLayerOnshore
@mixin WaterDepth, FaciesBase, SedimentBuffer, ActiveLayer
using ..Common
using ...Stencil: stencil

struct Facies <: AbstractFacies
    onshore_velocity
end

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
```

``` {.julia file=src/Models/OnshoreTransport.jl}
@compose module OnshoreTransport
@mixin Tag, H5Writer, CAProduction, ActiveLayerOnshore

using ..Common
using ..CAProduction: production
using ..TimeIntegration
using ..WaterDepth
using ModuleMixins: @for_each
using .H5Writer: run_model

export Input, Facies

function initial_state(input::Input)
    ca_state = CellularAutomaton.initial_state(input)
    for _ in 1:20
        CellularAutomaton.step!(input)(ca_state)
    end

    sediment_height = zeros(Height, input.box.grid_size...)
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...)

    return State(
        step=0, sediment_height=sediment_height,
        sediment_buffer=sediment_buffer,
        ca=ca_state.ca, ca_priority=ca_state.ca_priority)
end

function step!(input::Input)
    step_ca! = CellularAutomaton.step!(input)
    disintegrate! = disintegration(input)
    produce = production(input)
    transport = transportation(input)

    function (state::State)
        if mod(state.step, input.ca_interval) == 0
            step_ca!(state)
        end

        p = produce(state)
        d = disintegrate!(state)

        active_layer = p .+ d
        sediment = transport(state, active_layer)

        push_sediment!(state.sediment_buffer, sediment ./ input.depositional_resolution .|> NoUnits)
        state.sediment_height .+= sum(sediment; dims=1)[1,:,:]
        state.step += 1

        return H5Writer.DataFrame(
            production = p,
            disintegration = d,
            deposition = sediment)
    end
end

function write_header(fid, input::AbstractInput)
    @for_each(P -> P.write_header(fid, input), PARENTS)
end

end
```
