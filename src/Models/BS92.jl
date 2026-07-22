# ~/~ begin <<docs/src/bosscher-1992.md#src/Models/BS92.jl>>[init]
@compose module BS92
@mixin Tag, Diagnostics, Output, Production

using ..Common
using ..Production: uniform_production
using ..TimeIntegration
using ..WaterDepth
using ModuleMixins: @for_each
using ...Output: Frame

export Input, Facies

function initial_state(input::Input)
    bathymetry = initial_topography(input)
    return State(0, bathymetry)
end

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

function write_header(input::AbstractInput, output::AbstractOutput)
    @for_each(P -> P.write_header(input, output), PARENTS)
end

end
# ~/~ end
