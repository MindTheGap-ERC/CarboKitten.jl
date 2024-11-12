# ~/~ begin <<docs/src/components/time.md#src/Components/TimeIntegration.jl>>[init]
@compose module TimeIntegration
using ..Common
using HDF5
export time, n_writes, time_axis

@kwdef struct Input <: AbstractInput
    time::TimeProperties
end

@kwdef mutable struct State <: AbstractState
    step::Int
end

State(_::AbstractInput) = State(0)

time(input::AbstractInput, state::AbstractState) = input.time.t0 + state.step * input.time.Δt

write_times(input::AbstractInput) = write_times(input.time)
write_times(time::TimeProperties) = (0:n_writes(time)) .* (time.Δt * time.write_interval) .+ time.t0

time_axis(input::AbstractInput) = time_axis(input.time)
time_axis(time::TimeProperties) = (0:n_writes(time)) .* (time.Δt * time.write_interval) .+ time.t0

n_writes(input::AbstractInput) = n_writes(input.time)
n_writes(time::TimeProperties) = div(time.steps, time.write_interval)

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
# ~/~ end
