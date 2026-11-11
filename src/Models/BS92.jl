# ~/~ begin <<docs/src/bosscher-1992.md#src/Models/BS92.jl>>[init]
module BS92
import ModuleMixins.Procedural: requires
import ..Tag, ..Diagnostics, ..Output, ..Production

using ..Common
using ..Production: uniform_production
using ..TimeIntegration
using ..WaterDepth
using ...Output: Frame

# The model composes the records of its components. These are defined in the
# model module itself, so `requires` is defined first, instead of using
# `@component`, which defines it at the end of the module.
requires(::Val{BS92}) = (Tag, Diagnostics, Output, Production)
const Input = record_type(:Input, BS92)
const State = record_type(:State, BS92)
const Facies = record_type(:Facies, BS92)

export Input, Facies

initial_state(input::Input) = Common.initial_state(input)

function initial_frame(input::Input)
    return Frame(production=zeros(Sediment,length(input.facies), input.box.grid_size...),
                  disintegration=zeros(Sediment,length(input.facies), input.box.grid_size...),
                  deposition=zeros(Sediment,length(input.facies), input.box.grid_size...))
end

function step!(input::Input)
    τ = uniform_production(input)
    subside! = subsider(input)

    function (state::State)
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
