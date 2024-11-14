# ~/~ begin <<docs/src/components/cellular-automata.md#src/Components/CellularAutomaton.jl>>[init]
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
# ~/~ end
