# ~/~ begin <<docs/src/components/cellular-automata.md#src/Components/CellularAutomaton.jl>>[init]
@compose module CellularAutomaton
    @mixin Boxes, FaciesBase
    using ..Common
	using ...BoundaryTrait: get_bounded
	using Base: IndexCartesian
	using Random


    # ~/~ begin <<docs/src/components/cellular-automata.md#ca-input>>[init]
    @kwdef struct Facies <: AbstractFacies
        viability_range::Tuple{Int,Int} = (4, 10)
        activation_range::Tuple{Int,Int} = (6, 10)
        active::Bool = true
		neighbour_radius::Int = 2
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



    function rules(BT, facies, ca_priority, ca, i::CartesianIndex)
    cell_facies = get_bounded(BT, ca, i)

    function neighbour_count(f, r)
        c = 0
        for dx in -r:r, dy in -r:r
            (dx == 0 && dy == 0) && continue
            j = i + CartesianIndex(dx, dy)
            c += (get_bounded(BT, ca, j) == f)
        end
        return c
    end

    if cell_facies == 0
        for f in ca_priority
            r = facies[f].neighbour_radius
            n = neighbour_count(f, r)
            (a, b) = facies[f].activation_range
            if a <= n && n <= b
                return f
            end
        end
        return 0
    else
        r = facies[cell_facies].neighbour_radius
        n = neighbour_count(cell_facies, r)
        (a, b) = facies[cell_facies].viability_range
        return (a <= n && n <= b) ? cell_facies : 0
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
        for i in eachindex(IndexCartesian(), state.ca)
            tmp[i] = rules(BT, facies_, p, state.ca, i)
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
