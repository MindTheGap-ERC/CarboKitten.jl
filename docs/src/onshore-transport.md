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
using StaticArrays: Size
using ...Stencil: stencil!
using Unitful: Length

export otransportation

struct Facies <: AbstractFacies
  onshore_velocity
end

function onshore_transport_stencil(box::Box{BT}, Δt, ν, sf::F, out, w::AbstractArray{T}, C::AbstractArray{U}) where {BT<:Boundary{2},F,T<:Length,U<:Length}
  Δx = box.phys_scale
  d = ν * Δt

  stencil!(BT, Size(3, 3), out, w, C) do w, C
    sv, ss = sf(w[2, 2])

    adv = - (w[3, 2] - w[1, 2]) / (2Δx) * (d * (C[3, 2] - C[1, 2]) / (2Δx) + C[2, 2] * ss[1] * Δt) -
            (w[2, 3] - w[2, 1]) / (2Δx) * (d * (C[2, 3] - C[2, 1]) / (2Δx) + C[2, 2] * ss[2] * Δt)

    dif = - d * C[2, 2] * (w[3, 2] + w[2, 3] + w[1, 2] +
                           w[2, 1] - 4 * w[2, 2]) / (Δx)^2

    prd = - sv[1] * (C[3, 2] - C[1, 2]) * Δt / (2Δx) - sv[2] * (C[2, 3] - C[2, 1]) * Δt / (2Δx) + C[2, 2]

    return max(0.0u"m", adv + dif + prd)
  end
end

function odisintegration(input)
  max_h = input.disintegration_rate * input.time.Δt
  w = water_depth(input)
  output = Array{Float64,3}(undef, n_facies(input), input.box.grid_size...)

  return function (state)
    wn = w(state)
    h = min.(max_h, state.sediment_height)
    h[wn.<=0.0u"m"] .= 0.0u"m"
    state.sediment_height .-= h
    pop_sediment!(state.sediment_buffer, h ./ input.depositional_resolution .|> NoUnits, output)
    return output .* input.depositional_resolution
  end
end

function otransportation(input)
  w = water_depth(input)

  # We always return this array
  transported_output = Array{Amount,3}(undef, n_facies(input), input.box.grid_size...)
  box = input.box
  Δt = input.time.Δt
  fs = input.facies

  return function (state, active_layer::Array{Amount,3})
    wd = w(state)

    for (i, f) in pairs(fs)
      onshore_transport_stencil(
        box, Δt, f.diffusion_coefficient, f.onshore_velocity,
        view(transported_output, i, :, :),
        wd, view(active_layer, i, :, :))
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
    disintegrate! = ActiveLayerOnshore.odisintegration(input)
    produce = production(input)
    transport = ActiveLayerOnshore.otransportation(input)

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

## Xi & Burgess 2022

Xi & Burgess use the following equation for the phase velocity of waves as a function of depth:

$$v(w) = \sqrt{\frac{\lambda g}k} {\rm tanh} (k w),$$

where $k = 2\pi/\lambda$ and $g$ is the gravitational acceleration. It could be that the square root should be extended over the ${\rm tanh}$ function, following the derivation on **Airy wave theory** on Wikipedia. It should be noted that this velocity is the phase-velocity of surface waves, given the total depth of the water. To evaluate the transport velocity at deeper levels, we need to look at the **Stokes drift**. This will effectively multiply the phase velocity with a factor $\exp(-kw)$. We'll leave the proportionality as an input parameter.

$$v(w) \propto \tanh{k w} \exp{-kw}$$


```julia
using GLMakie
using Unitful

const g = 9.8u"m/s^2"

v(A, k, w) = A * tanh(k * w) * exp(-k * w)

let
    fig = Figure()
    ax = Axis(fig[1,1], yreversed=true, xlabel="v [m/s]", ylabel="wd [m]")

    wd = LinRange(0, 50, 10000)u"m"
    lines!(ax, v.(1.0u"m/s" * 3.331, 0.05u"1/m", wd) / u"m/s" .|> NoUnits, wd / u"m" .|> NoUnits)

    fig
end
```
