# ~/~ begin <<docs/src/models/ca-with-production.md#src/Models/CAP.jl>>[init]
@compose module CAP
@mixin Tag, Diagnostics, Output, CAProduction

using ..Common
using ..CAProduction: production
using ..TimeIntegration
using ..WaterDepth
using ...Output: Frame
using ModuleMixins: @for_each

export Input, Facies

function initial_state(input::Input)
    ca_state = CellularAutomaton.initial_state(input)
    for _ in 1:20
        CellularAutomaton.step!(input)(ca_state)
    end

    state = State(
        step=0, bathymetry=initial_topography(input),
        ca=ca_state.ca, ca_priority=ca_state.ca_priority)

    return state
end

function initial_frame(input::Input)
    return Frame(production=zeros(Sediment,n_facies(input), input.box.grid_size...),
                  deposition=zeros(Sediment,n_facies(input), input.box.grid_size...))
end

function step!(input::Input)
    τ = production(input)
    step_ca = CellularAutomaton.step!(input)
    subside! = subsider(input)

    function (state::State)
        if mod(state.step, input.ca_interval) == 0
            step_ca(state)
        end

        prod = τ(state)
        Δη = sum(prod; dims=1)[1, :, :]
        state.bathymetry .+= Δη
        subside!(state)
        state.step += 1

        return Frame(
            production=prod,
            deposition=prod)
    end
end

function write_header(input::AbstractInput, output::AbstractOutput)
    @for_each(P -> P.write_header(input, output), PARENTS)
end
end
# ~/~ end
