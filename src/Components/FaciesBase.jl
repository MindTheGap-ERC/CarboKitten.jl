# ~/~ begin <<docs/src/components/facies.md#src/Components/FaciesBase.jl>>[init]
@component module FaciesBase
using ..Common
import ..Common: write_header

using HDF5
export n_facies

@kwdef struct Facies <: FaciesPart
    name::Union{String,Nothing} = nothing
end

@kwdef struct Input <: InputPart
    facies::Vector{Facies} = []
end

n_facies(input::AbstractInput) = length(input.facies)

function write_header(::Input, input::AbstractInput, output::AbstractOutput)
    set_attribute(output, "n_facies", n_facies(input))

    for (i, f) in enumerate(input.facies)
        if f.name !== nothing
            set_attribute(output, "facies$(i)/name", f.name)
        end
    end
end
end
# ~/~ end
