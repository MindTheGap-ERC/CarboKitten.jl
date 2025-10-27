# ~/~ begin <<docs/src/finite-difference-transport.md#examples/transport/test_model.jl>>[init]
module TestModel

using CarboKitten
using CarboKitten: AbstractInput, AbstractState
using CarboKitten.Transport.Advection: transport
using Unitful

@kwdef struct Input <: AbstractInput
    box::Box
    time::TimeProperties
    initial_state::Array{Float64}
    topography::Array{typeof(1.0u"m")}

    diffusivity = 0.0u"m/Myr"
    wave_velocity = _ -> ((0.0u"m/Myr", 0.0u"m/Myr"), (0.0u"1/Myr", 0.0u"1/Myr"))

    solver

    diagnostics::Bool = false
end

@kwdef mutable struct State <: AbstractState
    time::typeof(1.0u"Myr")
    value::Array{Float64}
end

initial_state(input) = State(
    time = input.time.t0,
    value = copy(input.initial_state))

function step!(input)
    function (state)
        input.solver(
            (a, t) -> transport(
                input.box, input.diffusivity, input.wave_velocity,
                a, .-input.topography),
            state.value, state.time, input.time.Δt)
        state.time += input.time.Δt
    end
end

end
# ~/~ end
