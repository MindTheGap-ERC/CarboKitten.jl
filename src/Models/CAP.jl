# ~/~ begin <<docs/src/models/ca-with-production.md#src/Models/CAP.jl>>[init]
module CAP
import ModuleMixins.Procedural: requires
import ..Tag, ..Diagnostics, ..Output, ..CAProduction

using ..Common
using ..FaciesBase: n_facies
import ..CellularAutomaton
using ..CAProduction: production
using ..TimeIntegration
using ..WaterDepth
using ...Output: Frame

# The model composes the records of its components. These are defined in the
# model module itself, so `requires` is defined first, instead of using
# `@component`, which defines it at the end of the module.
requires(::Val{CAP}) = (Tag, Diagnostics, Output, CAProduction)
const Input = record_type(:Input, CAP)
const State = record_type(:State, CAP)
const Facies = record_type(:Facies, CAP)

export Input, Facies

initial_state(input::Input) = Common.initial_state(input)

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
end
# ~/~ end
