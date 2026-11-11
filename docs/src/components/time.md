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
    using CarboKitten.Components: TimeIntegration
    using CarboKitten.Components.TimeIntegration: write_times, time, n_writes
    Input = record_type(:Input, TimeIntegration)
    State = record_type(:State, TimeIntegration)

    let input = Input(time=TimeProperties(
        Δt = 0.2u"Myr", steps = 5))
      @test write_times(input) |> collect == [0.0, 0.2, 0.4, 0.6, 0.8, 1.0] .* u"Myr"
      @test time(input, State(step=4)) == 0.8u"Myr"
      @test n_writes(input) == 5
    end

    let input = Input(time=TimeProperties(
        Δt = 0.02u"Myr", steps = 50, t0 = -1.0u"Myr"))
      @test n_writes(input) == 50
      @test time(input, State(step=25)) == -0.5u"Myr"
    end
end
end
```

``` {.julia file=src/Components/TimeIntegration.jl}
@component module TimeIntegration
using ..Common
import ..Common: write_header
import ...CarboKitten: time_axis, n_writes, n_steps

using HDF5
export time, n_writes, time_axis

@kwdef struct Input <: InputPart
    time::TimeProperties
end

@kwdef mutable struct State <: StatePart
    step::Int
end

State(_::AbstractInput) = State(0)

Common.initial_state(::Input, input::AbstractInput) = State(input)

"""
    time(input, state)

Return the time given input and state.
"""
time(input::AbstractInput, state::AbstractState) = input.time.t0 + state.step * input.time.Δt
time(tprop::TimeProperties, state::AbstractState) = tprop.t0 + state.step * tprop.Δt
time(input::AbstractInput) =
    let t0 = input.time.t0, dt = input.time.Δt
        s -> t0 + s.step * dt
    end

write_times(input::AbstractInput) = write_times(input.time)
write_times(time::TimeProperties) = (0:n_writes(time)) .* time.Δt .+ time.t0

time_axis(input::AbstractInput) = time_axis(input.time)
n_writes(input::AbstractInput) = n_writes(input.time)
n_steps(input::AbstractInput) = input.time.steps

function write_header(::Input, input::AbstractInput, output::AbstractOutput)
    set_attribute(output, "t", write_times(input) .|> in_units_of(u"Myr"))
    set_attribute(output, "t0", input.time.t0 |> in_units_of(u"Myr"))
    set_attribute(output, "delta_t", input.time.Δt |> in_units_of(u"Myr"))
    set_attribute(output, "time_steps", input.time.steps)
end

end
```
