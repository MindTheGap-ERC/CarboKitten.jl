# ~/~ begin <<docs/src/onshore-transport.md#src/Components/ActiveLayerOnshore.jl>>[init]
@compose module ActiveLayerOnshore
@mixin WaterDepth, FaciesBase, SedimentBuffer, ActiveLayer
using ..Common
using ...Stencil: stencil
export otransportation

struct Facies <: AbstractFacies
  onshore_velocity
end

function onshore_transport_stencil(box::Box{BT}, Δt, ν, sf::F) where {BT<:Boundary{2},F}
  Δx = box.phys_scale

  function kernel(x)
    (w, P) = x[2, 2]
    sv, ss = sf(w)
    d = ν * Δt

    adv = (x[3, 2][1] - x[1, 2][1]) / (2Δx) * (d * (x[3, 2][2] - x[1, 2][2]) / (2Δx) - P * ss[1] * Δt) +
          (x[2, 3][1] - x[2, 1][1]) / (2Δx) * (d * (x[2, 3][2] - x[2, 1][2]) / (2Δx) - P * ss[2] * Δt)

    dif = d * P * (x[3, 2][1] + x[2, 3][1] + x[1, 2][1] +
                   x[2, 1][1] - 4 * x[2, 2][1]) / (Δx)^2

    prd = - sv[1] * (x[3, 2][2] - x[1, 2][2]) * Δt / (2Δx) - sv[2] * (x[2, 3][2] - x[2, 1][2]) * Δt / (2Δx) + x[2, 2][2]

    return max(0.0u"m", adv + dif + prd)
  end

  stencil(Tuple{Amount,Amount}, Amount, BT, (3, 3), kernel)
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
  stencils = [
    let stc = onshore_transport_stencil(input.box, input.time.Δt, f.diffusion_coefficient, f.onshore_velocity)
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
