# Time

``` {.julia file=src/Components/TimeIntegration.jl}
@compose module TimeIntegration
using ..Common
using HDF5
export time, n_writes

@kwdef struct Input <: AbstractInput
    time::TimeProperties
end

@kwdef mutable struct State <: AbstractState
    step::Int
end

State(input::AbstractInput) = State(0)

time(input::AbstractInput, state::AbstractState) = input.time.t0 + state.step * input.time.Δt

write_times(input::AbstractInput) = (0:n_writes(input)) .* (input.time.Δt * input.write_interval) .+ input.time.t0

n_writes(input::AbstractInput) = div(input.time.steps, input.time.write_interval)

function write_header(fid, input::AbstractInput)
    gid = fid["input"]
    attr = attributes(gid)

    gid["t"] = write_times(input) .|> in_units_of(u"Myr")
    attr["t0"] = input.time.t0 |> in_units_of(u"Myr")
    attr["delta_t"] = input.time.Δt |> in_units_of(u"Myr")
    attr["write_interval"] = input.time.write_interval
    attr["time_steps"] = input.time.steps
end

end
```
