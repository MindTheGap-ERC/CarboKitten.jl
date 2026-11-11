# ~/~ begin <<docs/src/models/alcap.md#src/Models/ALCAP.jl>>[init]
module ALCAP
import ModuleMixins.Procedural: requires
import ..Tag, ..Output, ..CAProduction, ..CAFeedback, ..ActiveLayer, ..InitialSediment, ..Diagnostics

using ..Common
using ..FaciesBase: n_facies
using ..SedimentBuffer: push_sediment, pop_sediment, initial_sediment_buffer
using ..ActiveLayer: lithification_factor
import ..CellularAutomaton
using ..CAProduction: production
using ..TimeIntegration
using ..WaterDepth
using ...Output: Frame

# The model composes the records of its components. These are defined in the
# model module itself, so `requires` is defined first, instead of using
# `@component`, which defines it at the end of the module.
requires(::Val{ALCAP}) = (Tag, Output, CAProduction, CAFeedback, ActiveLayer, InitialSediment, Diagnostics)
const Input = record_type(:Input, ALCAP)
const State = record_type(:State, ALCAP)
const Facies = record_type(:Facies, ALCAP)

export Input, Facies, BenthicProduction, PelagicProduction

# The optional state type `S` allows the same model to be run with a different
# composition, for instance `record_type(:State, ALCAP, CircularSedimentBuffer)`.
function initial_state(input::AbstractInput, ::Type{S}=State) where {S <: AbstractState}
    ca_state = CellularAutomaton.initial_state(input)
    for _ in 1:20
        CellularAutomaton.step!(input)(ca_state)
    end

    bathymetry = initial_topography(input)
    sediment_buffer = initial_sediment_buffer(input)
    sediment_thickness = zeros(Height, input.box.grid_size...)
    active_layer = zeros(Amount, n_facies(input), input.box.grid_size...)

    state = S(
        step=0,
        bathymetry=bathymetry,
        sediment_thickness=sediment_thickness,
        sediment_buffer=sediment_buffer,
        active_layer=active_layer,
        ca=ca_state.ca, ca_priority=ca_state.ca_priority)

    InitialSediment.push_initial_sediment!(input, state)
    return state
end

function initial_frame(input::AbstractInput)
    dep = stack(InitialSediment.initial_sediment(input.box, f) for f in input.facies; dims=1)
    return Frame(production=zeros(Sediment,size(dep)),
                  disintegration=zeros(Sediment,size(dep)),
                  deposition=dep)
end

function step!(input::AbstractInput)
    step_ca! = CellularAutomaton.step!(input)
    subside! = subsider(input)
    disintegrate! = ActiveLayer.disintegrator(input)
    produce = production(input)
    kill_unproductivity! = CAFeedback.ca_feedback(input)
    transport! = ActiveLayer.transporter(input)
    local_water_depth = water_depth(input)
    na = [CartesianIndex()]
    pf = lithification_factor(input)
    dtf = input.disintegration_transfer
    debug = input.diagnostics
    push! = push_sediment(input)
    pop! = pop_sediment(input)

    function (state::AbstractState)
        if debug
            @debug "step: " state.step
            @debug "   current active layer ambitus: " extrema(state.active_layer)
        end

        if mod(state.step, input.ca_interval) == 0
            step_ca!(state)
        end

        wd = local_water_depth(state)
        p = produce(state, wd)
        kill_unproductivity!(state.ca, p)
        d = disintegrate!(state)

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
        push!(state, deposit)
        state.active_layer .-= deposit
        subside!(state)
        state.step += 1

        return Frame(
            production=p,
            disintegration=d,
            deposition=deposit)
    end
end

include("ALCAP/Example.jl")

end
# ~/~ end
