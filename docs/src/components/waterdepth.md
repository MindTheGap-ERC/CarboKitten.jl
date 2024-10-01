
``` {.julia file=src/Components/WaterDepth.jl}
@compose module WaterDepth
    @mixin TimeIntegration, Boxes
    using ..Common

    export water_depth

    @kwdef struct Input <: AbstractInput
        sea_level          # function (t::Time) -> Length
        bedrock_elevation  # function (x::Location, y::Location) -> Length
        subsidence_rate::Rate
    end

    mutable struct State <: AbstractState
        sediment_height::Matrix{Height}
    end

    State(input::AbstractInput) = State(zeros(Height, input.box.grid_size...))

    function water_depth(input::AbstractInput)
        x, y = axes(input.box)
        eta0 = input.bedrock_elevation.(x, y')

        return function (state::AbstractState)
            return input.sea_level(state.time) .- eta0 .+
                (input.subsidence_rate * state.time) .- state.sediment_height
        end
    end
end
```
