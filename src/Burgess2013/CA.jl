# ~/~ begin <<docs/src/carbocat-ca.md#src/Burgess2013/CA.jl>>[init]
module CA

using ...BoundaryTrait
using ...Stencil
using ...Config: Box

export run_ca

# ~/~ begin <<docs/src/carbocat-ca.md#cycle-permutation>>[init]
cycle_permutation(n_species::Int) =
    (circshift(1:n_species, x) for x in Iterators.countfrom(0))
# ~/~ end
# ~/~ begin <<docs/src/carbocat-ca.md#burgess2013-rules>>[init]
function rules(facies::Vector{F}) where {F}
    function (neighbourhood::Matrix{Int}, order::Vector{Int})
        cell_facies = neighbourhood[3, 3]
        neighbour_count(f) = sum(neighbourhood .== f)
        if cell_facies == 0
            for f in order
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
end

# ~/~ begin <<docs/src/carbocat-ca.md#ca-stateful>>[init]
function step_ca(box::Box{BT}, facies) where {BT<:Boundary{2}}
    """Creates a propagator for the state, updating the celullar automaton in place.

    Contract: the `state` should have `ca::Matrix{Int}` and `ca_priority::Vector{Int}`
    members."""
    r = rules(facies)
    tmp = Matrix{Int}(undef, box.grid_size)
    stencil_op = stencil(Int, BT, (5, 5), r)

    function (state)
        stencil_op(state.ca, tmp, state.ca_priority)
        state.ca, tmp = tmp, state.ca
        state.ca_priority = circshift(state.ca_priority, 1)
        return state
    end
end
# ~/~ end

function run_ca(::Type{B}, facies::Vector{F}, init::Matrix{Int}, n_species::Int) where {B<:Boundary{2},F}
    r = rules(facies)
    Channel{Matrix{Int}}() do ch
        target = Matrix{Int}(undef, size(init))
        put!(ch, init)
        stencil_op = stencil(Int, B, (5, 5), r)
        for perm in cycle_permutation(n_species)
            stencil_op(init, target, perm)
            init, target = target, init
            put!(ch, init)
        end
    end
end
# ~/~ end

end
# ~/~ end