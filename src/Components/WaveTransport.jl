# ~/~ begin <<docs/src/onshore-transport.md#src/Components/WaveTransport.jl>>[init]
@compose module WaveTransport
@mixin WaterDepth, FaciesBase, SedimentBuffer
using ..Common
using ...Transport.Advection: transport!

struct Input <: AbstractInput
    disintegration_rate::typeof(1.0u"m/Myr")
    transport_solver = nothing
    transport_substeps::Int = 1
end

struct Facies <: AbstractFacies
    diffusivity
    wave_velocity
end

function disintegrator(input)
    max_h = input.disintegration_rate * input.time.Δt
    w = water_depth(input)
    output = Array{Float64,3}(undef, n_facies(input), input.box.grid_size...)
    depositional_resolution = input.depositional_resolution

    return function (state)
        wn = w(state)
        h = min.(max_h, state.sediment_height)
        h[wn.<=0.0u"m"] .= 0.0u"m"
        state.sediment_height .-= h
        pop_sediment!(state.sediment_buffer, h ./ depositional_resolution .|> NoUnits, output)
        return output .* depositional_resolution
    end
end

function transporter(input)
    solver = if input.transport_solver === nothing
        runge_kutta_4(input.box)
    else
        input.transport_solver
    end

    w = water_depth(input)
    box = input.box
    Δt = input.time.Δt / input.transport_substeps
    steps = input.transport_substeps
    fs = input.facies

    return function (state, C::Array{Amount,3})
        wd = w(state)

        for (i, f) in pairs(fs)
            for j in 1:steps
                solver(
                    (C, _) -> transport(
                        input.box, f.diffusivity, f.wave_velocity,
                        view(active_layer, i, :, :), wd),
                    C, state.time, Δt)
            end
        end
    end
end

end
# ~/~ end
