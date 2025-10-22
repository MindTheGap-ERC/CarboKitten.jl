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
    sediment_height = zeros(Height, input.box.grid_size...)
    return State(0, sediment_height)
end

function step!(input::Input)
    τ = uniform_production(input)
    function (state::State)
        prod = τ(state)
        Δη = sum(prod; dims=1)[1, :, :]
        state.sediment_height .+= Δη
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
