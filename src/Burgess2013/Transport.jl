# ~/~ begin <<docs/src/carbocat-transport.md#src/Burgess2013/Transport.jl>>[init]
module Transport

using ...BoundaryTrait
using ...Stencil
using Transducers

struct Product
    species::Int
    amount::Float64
end

Base.zero(::Type{Product}) = Product(0, 0.0)

struct Deposit{N}
    amount::NTuple{N, Float64}
end

Base.zero(::Type{Deposit{N}}) where {N} =
    Deposit{N}(ntuple(_ -> zero(Float64), N))

function deposit(
        ::Type{B},
        production::Matrix{Product},
        elevation::Matrix{Float64},
        transported::Matrix{Deposit{N}},
        lim::Float64,
        idx::CartesianIndex,
        p::Product
    ) where {B <: Boundary{2}, N}

    if p.amount <= lim
        transported[idx][p.species] += p.amount
        return
    end

    shape = size(production)
    targets = CartesianIndices((-1:1,-1:1)) |>
        Filter(Δi -> Δi!=CartesianIndex(0,0)) |>
        Map(Δi -> offset_index(B, shape, idx, Δi)) |>
        Filter(j -> !isnothing(j) &&
                    elevation[j] >= elevation[idx] &&
                    production[j].species == 0) |>
        collect

    if isempty(targets)
        transported[idx][p.species] += p.amount
        return
    end
    transported[idx][p.species] += p.amount / 2
    for j in targets
        q = Product(p.species, p.amount / (2 * length(targets)))
        deposit(B, production, elevation, transported, lim, j, q)
    end
end

function transport(
        ::Type{B},
        n_species::Int,
        production::Matrix{Product},
        elevation::Matrix{Float64},
        fraction::Float64,
        lim::Float64
    ) where {B <: Boundary{2}}

    shape = size(production)
    result = zeros(Transported{n_species})

    for i in CartesianIndices(shape)
        production[i].amount <= lim && continue
        p = copy(production[i])
        p.amount *= fraction
        production[i].amount -= p.amount
        deposit(B, production, elevation, result, lim, i, p)
    end

    return result
end

end
# ~/~ end