# ~/~ begin <<docs/src/bosscher-1992.md#src/Model/BS92.jl>>[init]
@compose module BS92
@mixin H5Writer, Production

using ..Common
using ..Production: uniform_production
using ..TimeIntegration
using ..WaterDepth

export Input, Facies, run

function initial_state(input::Input)
    sediment_height = zeros(Height, input.box.grid_size...)
    return State(0, sediment_height)
end

function step!(input::Input)
    τ = uniform_production(input)
    function (state::State)
        prod = τ(state) .* input.time.Δt
        Δη = sum(prod; dims=1)[1, :, :]
        state.sediment_height .+= Δη
        state.step += 1
        return H5Writer.DataFrame(
            production = prod,
            deposition = prod)
    end
end

function write_header(fid, input::AbstractInput)
    @for_each(P -> P.write_header(fid, input), PARENTS)
end

end
# ~/~ end
