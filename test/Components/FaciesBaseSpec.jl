# ~/~ begin <<docs/src/components/facies.md#test/Components/FaciesBaseSpec.jl>>[init]
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
# ~/~ end
