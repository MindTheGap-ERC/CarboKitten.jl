# Time

``` {.julia file=src/Components/TimeIntegration.jl}
@compose module TimeIntegration
using ..Common
using HDF5

@kwdef struct Input <: AbstractInput
    time::TimeProperties
end

mutable struct State <: AbstractState
    step::Int
end

State(input::AbstractInput) = State(0)

time(input::AbstractInput, state::AbstractState) = state.step * input.time.Δt

function write_header(fid, input::AbstractInput)
    t = (0:input.time.steps) .* input.time.Δt
    gid = fid["input"]
    attr = attributes(gid)

    gid["t"] = t .|> in_units_of(u"Myr")
    attr["delta_t"] = input.time.Δt |> in_units_of(u"Myr")
    attr["write_interval"] = input.time.write_interval
    attr["time_steps"] = input.time.steps
end
end
```
