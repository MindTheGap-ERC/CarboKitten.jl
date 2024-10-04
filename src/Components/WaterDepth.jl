# ~/~ begin <<docs/src/components/waterdepth.md#src/Components/WaterDepth.jl>>[init]
@compose module WaterDepth
    @mixin TimeIntegration, Boxes
    using ..Common
    using CarboKitten.Config: axes

    export water_depth

    @kwdef struct Input <: AbstractInput
        sea_level          # function (t::Time) -> Length
        bedrock_elevation  # function (x::Location, y::Location) -> Length
        subsidence_rate::Rate
    end

    mutable struct State <: AbstractState
        sediment_height::Matrix{Height}
    end

    function water_depth(input::AbstractInput)
        x, y = axes(input.box)
        eta0 = input.bedrock_elevation.(x, y')

        return function (state::AbstractState)
            t = TimeIntegration.time(input, state)
            return input.sea_level(t) .- eta0 .+
                (input.subsidence_rate * t) .- state.sediment_height
        end
    end
end
# ~/~ end