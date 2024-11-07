# Cellular Automata

```component-dag
CarboKitten.Components.CellularAutomaton
```

``` {.julia file=test/Components/CellularAutomatonSpec.jl}
module CellularAutomatonSpec

using Test
using CarboKitten.Components.Common
using CarboKitten.Components: CellularAutomaton as CA

@testset "Components/CellularAutomaton" begin
    let facies = fill(CA.Facies((4, 10), (6, 10)), 3),
        input1 = CA.Input(
            box=Box{Periodic{2}}(grid_size=(50, 50), phys_scale=1.0u"m"),
            facies=facies),
        input2 = CA.Input(
            box=Box{Periodic{2}}(grid_size=(50, 50), phys_scale=1.0u"m"),
            facies=facies,
            ca_random_seed=1)

        state1 = CA.initial_state(input1)
        state2 = CA.initial_state(input2)
        state3 = CA.initial_state(input2)

        @test CA.initial_state(input1).ca == CA.initial_state(input1).ca
        @test state1.ca != state2.ca

        step! = CA.step!(input1)  # inputs have same rules
        for i in 1:20
            step!(state1)
            step!(state2)
            step!(state3)
        end

        @test state1.ca != state2.ca
        @test state2.ca == state3.ca
    end
end

end
```

``` {.julia file=src/Components/CellularAutomaton.jl}
@compose module CellularAutomaton
    @mixin Boxes, FaciesBase
    using ..Common
    using Random
    using ...Burgess2013.CA: step_ca

    @kwdef struct Facies <: AbstractFacies
        viability_range::Tuple{Int,Int} = (4, 10)
        activation_range::Tuple{Int,Int} = (6, 10)
    end

    @kwdef struct Input <: AbstractInput
        ca_interval::Int      = 1
        ca_random_seed::Int   = 0
    end

    @kwdef mutable struct State <: AbstractState
        ca::Matrix{Int}
        ca_priority::Vector{Int}
    end

    function initial_state(input::AbstractInput)
        n_facies = length(input.facies)
        ca = rand(MersenneTwister(input.ca_random_seed), 0:n_facies, input.box.grid_size...)
        return State(ca, 1:n_facies |> collect)
    end

    function step!(input::AbstractInput)
        return step_ca(input.box, input.facies)
    end
end
```
