# ~/~ begin <<docs/src/components/time.md#src/Components/TimeIntegration.jl>>[init]
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

time(input::AbstractInput, state::AbstractState) = state.step * input.time.Δt

n_writes(input::AbstractInput) = div(input.time.steps, input.time.write_interval)

function write_header(fid, input::AbstractInput)
    t = (0:n_writes(input)) .* (input.time.Δt * input.time.write_interval)
    gid = fid["input"]
    attr = attributes(gid)

    gid["t"] = t .|> in_units_of(u"Myr")
    attr["delta_t"] = input.time.Δt |> in_units_of(u"Myr")
    attr["write_interval"] = input.time.write_interval
    attr["time_steps"] = input.time.steps
end

end
# ~/~ end