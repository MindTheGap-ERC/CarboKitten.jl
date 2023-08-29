# ~/~ begin <<docs/src/carbocat-ca.md#src/Burgess2013/CA.jl>>[init]
module CA

using MindTheGap.Stencil

# ~/~ begin <<docs/src/carbocat-ca.md#cycle-permutation>>[init]
cycle_permutation(n_species::Int) =
    (circshift(1:n_species, x) for x in Iterators.countfrom(0))
# ~/~ end
# ~/~ begin <<docs/src/carbocat-ca.md#burgess2013-rules>>[init]
function rules(neighbourhood::Matrix{Int}, order::Vector{Int})
    cell_species = neighbourhood[3, 3]
    neighbour_count(species) = sum(neighbourhood .== species)
    if cell_species == 0
        for species in order
            n = neighbour_count(species)
            if 6 <= n && n <= 10
                return species
            end
        end
        0
    else
        n = neighbour_count(cell_species) - 1
        (4 <= n && n <= 10 ? cell_species : 0)
    end
end
# ~/~ end

function run(::Type{B}, init::Matrix{Int}, n_species::Int) where {B <: Boundary{2}}
    Channel{Matrix{Int}}() do ch
        target = Matrix{Int}(undef, size(init))
        put!(ch, init)
        stencil_op = stencil(Int, B, (5, 5), rules)
        for perm in cycle_permutation(n_species)
            stencil_op(init, target, perm)
            init, target = target, init
            put!(ch, init)
        end
    end
end

end
# ~/~ end