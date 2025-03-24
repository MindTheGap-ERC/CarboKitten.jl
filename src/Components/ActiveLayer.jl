# ~/~ begin <<docs/src/active-layer-transport.md#src/Components/ActiveLayer.jl>>[init]
@compose module ActiveLayer
@mixin WaterDepth, FaciesBase, SedimentBuffer

export disintegrator, transporter

using ..Common
using CarboKitten.Transport.Advection: transport
using CarboKitten.Transport.Solvers: runge_kutta_4
using Unitful

@kwdef struct Facies <: AbstractFacies
    diffusion_coefficient::typeof(1.0u"m/yr") = 0.0u"m/Myr"
    wave_velocity = _ -> ((0.0u"m/Myr", 0.0u"m/Myr"), (0.0u"1/Myr", 0.0u"1/Myr"))
end

@kwdef struct Input <: AbstractInput
    disintegration_rate::Rate = 50.0u"m/Myr"
    transport_solver = nothing
    transport_substeps::Int = 1 
end

"""
    disintegrator(input) -> f!

Prepares the disintegration step. Returns a function `f!(state::State)`. The returned function
modifies the state, popping sediment from the `sediment_buffer` and returns an array of `Amount`.
"""
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

"""
    transporter(input::Input) -> f

Prepares the transportation step. Returns a function `f(state::State, active_layer)`,
transporting the active layer, returning a transported `Amount` of sediment.
"""
function transporter(input)
    solver = if input.transport_solver === nothing
        runge_kutta_4(typeof(1.0u"m"), input.box)
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
                        input.box, f.diffusion_coefficient, f.wave_velocity,
                        C, wd),
                    view(C, i, :, :), TimeIntegration.time(input, state), Δt)
            end
        end
    end
end

end
# ~/~ end
