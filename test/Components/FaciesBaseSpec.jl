# ~/~ begin <<docs/src/components/facies.md#test/Components/FaciesBaseSpec.jl>>[init]
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
# ~/~ end