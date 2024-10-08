# ~/~ begin <<docs/src/bosscher-1992.md#src/Model/BS92.jl>>[init]
@compose module BS92
@mixin Production

using ..Common
using CSV
using DataFrames
using Interpolations
using ..Production: uniform_production
using ..TimeIntegration
using ..WaterDepth

function State(input::Input)
    sediment_height = zeros(Height, input.box.grid_size...)
    return State(0, sediment_height)
end

function step(input::Input)
    τ = uniform_production(input)
    function (state::State)
        prod = τ(state) .* input.time.Δt
        Δη = sum(prod; dims=1)[1, :, :]
        state.sediment_height .+= Δη
        state.step += 1
        return prod
    end
end

function sealevel_curve()
    data = DataFrame(CSV.File("../data/bs92-sealevel-curve.csv"))
    linear_interpolation(data.time, data.depth)
end

@kwdef struct Frame
    deposition::Array{Amount,3}
end

function write_header(fid, input::AbstractInput)
    create_group(fid, "input")
    @for_each(P -> P.write_header(fid, input), PARENTS)
end

function run(input::Input)
    step! = step(input)
    state = State(input)

    n_writes = input.time.steps ÷ input.time.write_interval
    Channel{Frame}() do ch
        for i = 1:n_writes
            prod = zeros(Amount, n_facies(input), input.box.grid_size...)
            for _ = 1:input.time.write_interval
                prod .+= step!(state)
            end
            f = Frame(
                deposition=prod,
                sediment_height=copy(state.sediment_height))
            put!(ch, f)
        end
    end
end
end
# ~/~ end
