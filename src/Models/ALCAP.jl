# ~/~ begin <<docs/src/model-alcap.md#src/Models/ALCAP.jl>>[init]
@compose module ALCAP
@mixin Tag, Output, CAProduction, ActiveLayer, InitialSediment

using ..Common
using ..CAProduction: production
using ..TimeIntegration
using ..WaterDepth: water_depth
using ...Output: Frame
using ModuleMixins: @for_each

export Input, Facies

function initial_state(input::AbstractInput)
    ca_state = CellularAutomaton.initial_state(input)
    for _ in 1:20
        CellularAutomaton.step!(input)(ca_state)
    end

    sediment_height = zeros(Height, input.box.grid_size...)
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...)
    active_layer = zeros(Amount, n_facies(input), input.box.grid_size...)

    state = State(
        step=0, sediment_height=sediment_height,
        sediment_buffer=sediment_buffer,
        active_layer=active_layer,
        ca=ca_state.ca, ca_priority=ca_state.ca_priority)

    InitialSediment.push_initial_sediment!(input, state)
    return state
end

function step!(input::Input)
    step_ca! = CellularAutomaton.step!(input)
    disintegrate! = ActiveLayer.disintegrator(input)
    produce = production(input)
    transport! = ActiveLayer.transporter(input)
    local_water_depth = water_depth(input)
    na = [CartesianIndex()]
    pf = cementation_factor(input)
    dtf = input.disintegration_transfer

    function (state::State)
        if mod(state.step, input.ca_interval) == 0
            step_ca!(state)
        end

        wd = local_water_depth(state)
        p = produce(state, wd)
        d = disintegrate!(state)

        state.active_layer .+= p
        state.active_layer .+= dtf(d)
        transport!(state)

        deposit = pf .* state.active_layer
        push_sediment!(state.sediment_buffer, deposit ./ input.depositional_resolution .|> NoUnits)
        state.active_layer .-= deposit
        state.sediment_height .+= sum(deposit; dims=1)[1, :, :]
        state.step += 1

        return Frame(
            production=p,
            disintegration=d,
            deposition=deposit)
    end
end

function write_header(input::AbstractInput, output::AbstractOutput)
    @for_each(P -> P.write_header(input, output), PARENTS)
end

include("ALCAP/Example.jl")

end
# ~/~ end
