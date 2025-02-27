# ~/~ begin <<docs/src/onshore-transport.md#src/Components/ActiveLayerOnshore.jl>>[init]
@compose module ActiveLayerOnshore
@mixin WaterDepth, FaciesBase, SedimentBuffer, ActiveLayer
using ..Common
using StaticArrays: Size
using ...Stencil: stencil!
export otransportation

struct Facies <: AbstractFacies
  onshore_velocity
end

function onshore_transport_stencil(box::Box{BT}, Δt, ν, sf::F, out, w, C) where {BT<:Boundary{2},F}
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
# ~/~ end
