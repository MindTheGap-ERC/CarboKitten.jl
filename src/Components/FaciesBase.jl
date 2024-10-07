# ~/~ begin <<docs/src/components/facies.md#src/Components/FaciesBase.jl>>[init]
@compose module FaciesBase
    using ..Common
    export n_facies

    @kwdef struct Facies <: AbstractFacies
    end

    @kwdef struct Input <: AbstractInput
        facies::Vector{Facies} = []
    end

    n_facies(input::AbstractInput) = length(input.facies)
end
# ~/~ end
