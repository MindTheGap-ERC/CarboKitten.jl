# ~/~ begin <<docs/src/bosscher-1992.md#src/Model/BS92.jl>>[init]
@compose module BS92
    @mixin UniformProduction

    using ..Common
    using CSV
    using DataFrames
    using Interpolations
    using ..UniformProduction: uniform_production
    using ..TimeIntegration
    using ..WaterDepth

    function State(input::Input)
        ti_state = TimeIntegration.State(input)
        sediment_height = zeros(Height, input.box.grid_size...)
        return State(0, sediment_height)
    end

    function step(input::Input)
        τ = uniform_production(input)
        function (state::State)
            prod = τ(state) .* input.time.Δt
            Δη = sum(prod; dims=1)[1,:,:]
            state.sediment_height .+= Δη
            state.step += 1
            return prod
        end
    end

    function sealevel_curve()
        data = DataFrame(CSV.File("../data/bs92-sealevel-curve.csv"))
        linear_interpolation(data.time, data.depth)
    end

    struct Frame
        deposition::Array{Amount, 3}
        sediment_height::Array{Amount, 2}
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
                put!(ch, Frame(prod, copy(state.sediment_height)))
            end
        end
    end
end
# ~/~ end
