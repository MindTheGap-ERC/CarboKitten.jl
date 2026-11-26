# Facies Base

```component-dag
CarboKitten.Components.FaciesBase
```

Base module for including facies.

``` {.julia file=test/Components/FaciesBaseSpec.jl}
module FaciesBaseSpec
    using Test
    using CarboKitten.Components.Common
    using CarboKitten.Components: FaciesBase
    using CarboKitten.Components.FaciesBase: n_facies

    const Facies = record_type(:Facies, FaciesBase)
    const Input = record_type(:Input, FaciesBase)

    @testset "Components/FaciesBase" begin
        let input = Input(facies=fill(Facies(), 23))
            @test n_facies(input) == 23
        end
    end
end
```

``` {.julia file=src/Components/FaciesBase.jl}
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
```
