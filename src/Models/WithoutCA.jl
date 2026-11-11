# ~/~ begin <<docs/src/models/without-ca.md#src/Models/WithoutCA.jl>>[init]
module WithoutCA
import ModuleMixins.Procedural: requires
import ..Tag, ..Diagnostics, ..Output, ..Production, ..ActiveLayer, ..InitialSediment

using ..Common
using ..FaciesBase: n_facies
using ..SedimentBuffer: push_sediment
using ..ActiveLayer: lithification_factor
using ..Production: uniform_production
using ..TimeIntegration
using ..WaterDepth
using ...Output: Frame

# The model composes the records of its components. These are defined in the
# model module itself, so `requires` is defined first, instead of using
# `@component`, which defines it at the end of the module.
requires(::Val{WithoutCA}) = (Tag, Diagnostics, Output, Production, ActiveLayer, InitialSediment)
const Input = record_type(:Input, WithoutCA)
const State = record_type(:State, WithoutCA)
const Facies = record_type(:Facies, WithoutCA)

export Input, Facies

initial_state(input::Input) = Common.initial_state(State, input)

function initial_frame(input::Input)
    dep = stack(InitialSediment.initial_sediment(input.box, f) for f in input.facies; dims=1)
    return Frame(production=zeros(Sediment,size(dep)),
                  disintegration=zeros(Sediment,size(dep)),
                  deposition=dep)
end

function step!(input::Input)
    disintegrate! = ActiveLayer.disintegrator(input)
    transport! = ActiveLayer.transporter(input)
    produce = uniform_production(input)
    dt = input.time.Δt
    local_water_depth = water_depth(input)
    na = [CartesianIndex()]
    pf = lithification_factor(input)
    dtf = input.disintegration_transfer
    push! = push_sediment(input)
    subside! = subsider(input)

    function (state::State)
        wd = local_water_depth(state)
        p = produce(state, wd)
        d = disintegrate!(state)

        state.active_layer .+= p
        state.active_layer .+= dtf(d)
        transport!(state)

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

end
# ~/~ end
