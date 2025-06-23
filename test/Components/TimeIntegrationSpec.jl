# ~/~ begin <<docs/src/components/time.md#test/Components/TimeIntegrationSpec.jl>>[init]
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
        Δt = 0.02u"Myr", steps = 50, t0 = -1.0u"Myr"))
      @test n_writes(input) == 50
      @test time(input, State(25)) == -0.5u"Myr"
    end
end
end
# ~/~ end
