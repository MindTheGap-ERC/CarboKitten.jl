# ~/~ begin <<docs/src/models/envcap.md#src/Models/EnvCAP.jl>>[init]
include("EnvCAP/EnvMapping.jl")

@compose module EnvCAP
@mixin Tag, Output, CAProduction, CAFeedback, ActiveLayer, InitialSediment, Diagnostics

using ..Common
using ..CAProduction: production
using ..TimeIntegration
using ..WaterDepth: water_depth
using ...Output: Frame
using ModuleMixins: @for_each
using Random

import ..EnvMapping
using ..EnvMapping: apply_prior_bias!

export Input, Facies, BenthicProduction, PelagicProduction
export EnvMapping

# ~/~ begin <<docs/src/models/envcap.md#envcap-input>>[init]
@kwdef struct Input <: AbstractInput
    factory_prior::Union{Array{Float64,4},Nothing} = nothing
    ca_refinement::Float64                         = 0.0
    env_random_seed::Int                           = 1
end
# ~/~ end
# ~/~ begin <<docs/src/models/envcap.md#envcap-initial-state>>[init]
function initial_state(input::AbstractInput)
    ca_state = CellularAutomaton.initial_state(input)
    for _ in 1:20
        CellularAutomaton.step!(input)(ca_state)
    end

    sediment_height = zeros(Height, input.box.grid_size...)
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...)
    active_layer    = zeros(Amount, n_facies(input), input.box.grid_size...)

    state = State(
        step            = 0,
        sediment_height = sediment_height,
        sediment_buffer = sediment_buffer,
        active_layer    = active_layer,
        ca              = ca_state.ca,
        ca_priority     = ca_state.ca_priority)

    InitialSediment.push_initial_sediment!(input, state)
    return state
end
# ~/~ end

function initial_frame(input::Input)
    dep = stack(InitialSediment.initial_sediment(input.box, f) for f in input.facies; dims=1)
    return Frame(production      = zeros(Sediment, size(dep)),
                 disintegration  = zeros(Sediment, size(dep)),
                 deposition      = dep)
end

# ~/~ begin <<docs/src/models/envcap.md#envcap-step>>[init]
function step!(input::Input)
    step_ca!           = CellularAutomaton.step!(input)
    disintegrate!      = ActiveLayer.disintegrator(input)
    produce            = production(input)
    kill_unproductive! = CAFeedback.ca_feedback(input)
    transport!         = ActiveLayer.transporter(input)
    local_water_depth  = water_depth(input)
    pf                 = lithification_factor(input)
    dtf                = input.disintegration_transfer
    debug              = input.diagnostics

    prior     = input.factory_prior
    α         = input.ca_refinement
    @assert 0.0 <= α <= 1.0 "ca_refinement must be between 0 and 1"
    env_rng   = MersenneTwister(input.env_random_seed)
    has_prior = prior !== nothing && α > 0.0

    function (state::State)
        if debug
            @debug "envcap step: " state.step
        end

        if mod(state.step, input.ca_interval) == 0
            step_ca!(state)
            if has_prior
                apply_prior_bias!(
                    state.ca,
                    prior,
                    state.sediment_height,
                    input.depositional_resolution,
                    α,
                    env_rng,
                )
            end
        end

        wd = local_water_depth(state)
        p  = produce(state, wd)
        kill_unproductive!(state.ca, p)
        d  = disintegrate!(state)

        state.active_layer .+= p
        state.active_layer .+= dtf(d)

        if debug
            @debug "   post-production ambitus: " extrema(state.active_layer)
        end

        transport!(state)

        if debug
            @debug "   post-transport ambitus: " extrema(state.active_layer)
        end

        deposit = pf .* state.active_layer
        push_sediment!(state.sediment_buffer, deposit ./ input.depositional_resolution .|> NoUnits)
        state.active_layer  .-= deposit
        state.sediment_height .+= sum(deposit; dims=1)[1, :, :]
        state.step += 1

        return Frame(
            production     = p,
            disintegration = d,
            deposition     = deposit)
    end
end
# ~/~ end

function write_header(input::AbstractInput, output::AbstractOutput)
    @for_each(P -> P.write_header(input, output), PARENTS)
end

end
# ~/~ end
