# ~/~ begin <<docs/src/components/cellular-automata.md#src/Components/CellularAutomaton.jl>>[init]
@compose module CellularAutomaton
    @mixin Boxes, FaciesBase
    using ..Common
    using ...Stencil: stencil!
    using Random

    # ~/~ begin <<docs/src/components/cellular-automata.md#ca-input>>[init]
    @kwdef struct Facies <: AbstractFacies
        viability_range::Tuple{Int,Int} = (4, 10)
        activation_range::Tuple{Int,Int} = (6, 10)
        active::Bool = true
    end

    @kwdef struct Input <: AbstractInput
        ca_interval::Int      = 1
        ca_random_seed::Int   = 0
    end
    # ~/~ end
    # ~/~ begin <<docs/src/components/cellular-automata.md#ca-state>>[init]
    @kwdef mutable struct State <: AbstractState
        ca::Matrix{Int}
        ca_priority::Vector{Int}
    end
    # ~/~ end
    # ~/~ begin <<docs/src/components/cellular-automata.md#ca-step>>[init]
    function rules(facies, ca_priority, neighbourhood)
        cell_facies = neighbourhood[3, 3]
        neighbour_count(f) = sum(neighbourhood .== f)
        if cell_facies == 0
            for f in ca_priority
                n = neighbour_count(f)
                (a, b) = facies[f].activation_range
                if a <= n && n <= b
                    return f
                end
            end
            0
        else
            n = neighbour_count(cell_facies) - 1
            (a, b) = facies[cell_facies].viability_range
            (a <= n && n <= b ? cell_facies : 0)
        end
    end
    # ~/~ end
    # ~/~ begin <<docs/src/components/cellular-automata.md#ca-step>>[1]
    """
        step_ca(box, facies)

    Creates a propagator for the state, updating the celullar automaton in place.

    Contract: the `state` should have `ca::Matrix{Int}` and `ca_priority::Vector{Int}`
    members.
    """
    function step_ca(box::Box{BT}, facies) where {BT<:Boundary{2}}
        tmp = Matrix{Int}(undef, box.grid_size)
        facies_ = facies

        function (state)
            p = state.ca_priority
            stencil!(BT, Size(5, 5), tmp, state.ca) do nb
                rules(facies_, p, nb)
            end
            state.ca, tmp = tmp, state.ca
            state.ca_priority = circshift(state.ca_priority, 1)
            return state
        end
    end
    # ~/~ end

    function initial_state(input::AbstractInput)
        n_facies = length(input.facies)
        active_facies = 1:n_facies |> filter(f->input.facies[f].active==true)
        ca = rand(MersenneTwister(input.ca_random_seed), [0; active_facies], input.box.grid_size...)
        return State(ca, active_facies |> collect)
    end

    function step!(input::AbstractInput)
        return step_ca(input.box, input.facies)
    end
end
# ~/~ end
