# Time

```component-dag
CarboKitten.Components.TimeIntegration
```

```@docs
CarboKitten.Components.TimeIntegration.time
```

``` {.julia file=test/Components/TimeIntegrationSpec.jl}
module TimeIntegrationSpec
using Test
using CarboKitten.Components.Common

@testset "Components/TimeIntegration" begin
    using CarboKitten.Components.TimeIntegration: write_times, Input, State, time, n_writes

    let input = Input(time=TimeProperties(
        Δt = 0.2u"Myr", steps = 5))
      @test write_times(input) |> collect == [0.0, 0.2, 0.4, 0.6, 0.8, 1.0] .* u"Myr"
      @test time(input, State(4)) == 0.8u"Myr"
      @test n_writes(input) == 5
    end

    let input = Input(time=TimeProperties(
        Δt = 0.02u"Myr", steps = 50, write_interval = 10, t0 = -1.0u"Myr"))
      @test write_times(input) |> collect == [-1.0, -0.8, -0.6, -0.4, -0.2, 0.0] .* u"Myr"
      @test n_writes(input) == 5
      @test time(input, State(25)) == -0.5u"Myr"
    end
end
end
```

``` {.julia file=src/Components/TimeIntegration.jl}
@compose module TimeIntegration
using ..Common
import ...CarboKitten: time_axis, n_writes

using HDF5
export time, n_writes, time_axis

@kwdef struct Input <: AbstractInput
    time::TimeProperties
end

@kwdef mutable struct State <: AbstractState
    step::Int
end

State(_::AbstractInput) = State(0)

"""
    time(input, state)

Return the time given input and state.
"""
time(input::AbstractInput, state::AbstractState) = input.time.t0 + state.step * input.time.Δt

write_times(input::AbstractInput) = write_times(input.time)
write_times(time::TimeProperties) = (0:n_writes(time)) .* (time.Δt * time.write_interval) .+ time.t0

time_axis(input::AbstractInput) = time_axis(input.time)
n_writes(input::AbstractInput) = n_writes(input.time)

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
