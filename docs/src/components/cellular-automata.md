# Cellular Automata

``` {.julia file=src/Components/CellularAutomata.jl}
@compose module CellularAutomata
    @mixin Boxes, FaciesConcept
    using ..Common
    using Random

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
    end
end
```
