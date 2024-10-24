# Facies Base

Base module for including facies.

``` {.julia file=test/Components/FaciesBaseSpec.jl}
module FaciesBaseSpec
    using Test
    using CarboKitten.Components.Common
    using CarboKitten.Components.FaciesBase: Facies, Input, n_facies

    @testset "Components/FaciesBase" begin
        let input = Input(facies=fill(Facies(), 23))
            @test n_facies(input) == 23
        end
    end
end
```

``` {.julia file=src/Components/FaciesBase.jl}
@compose module FaciesBase
using ..Common
using HDF5
export n_facies

@kwdef struct Facies <: AbstractFacies
end

@kwdef struct Input <: AbstractInput
    facies::Vector{Facies} = []
end

n_facies(input::AbstractInput) = length(input.facies)

function write_header(fid, input::AbstractInput)
    attr = attributes(fid["input"])
    attr["n_facies"] = n_facies(input)
    for i in 1:n_facies(input)
        create_group(fid["input"], "facies$(i)")
    end
end
end
```
