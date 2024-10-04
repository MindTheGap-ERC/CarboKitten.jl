# Cellular Automata

``` {.julia file=src/Components/CellularAutomaton.jl}
@compose module CellularAutomaton
    @mixin Boxes, FaciesBase
    using ..Common
    using Random
    using ...Burgess2013.CA: step_ca

    @kwdef struct Facies <: AbstractFacies
        viability_range::Tuple{Int,Int}
        activation_range::Tuple{Int,Int}
    end

    @kwdef struct Input <: AbstractInput
        ca_interval::Int      = 1
        ca_random_seed::Int   = 0
    end

    mutable struct State <: AbstractState
        ca::Matrix{Int}
        ca_priority::Vector{Int}
    end

    function initial_state(input::AbstractInput)
        n_facies = length(input.facies)
        ca = rand(MersenneTwister(input.ca_random_seed), 0:n_facies, input.box.grid_size...)
        return State(ca, 1:n_facies |> collect)
    end

    function stepper(input::AbstractInput)
        return step_ca(input.box, input.facies)
    end
end
```
