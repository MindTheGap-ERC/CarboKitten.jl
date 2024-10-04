# Time

``` {.julia file=src/Components/TimeIntegration.jl}
@compose module TimeIntegration
    using ..Common

    @kwdef struct Input <: AbstractInput
        time::TimeProperties
    end

    mutable struct State <: AbstractState
        step::Int
    end

    State(input::AbstractInput) = State(0)

    time(input::AbstractInput, state::AbstractState) = state.step * input.time.Δt
end
```
