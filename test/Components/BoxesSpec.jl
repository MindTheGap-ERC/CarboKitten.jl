# ~/~ begin <<docs/src/components/boxes.md#test/Components/BoxesSpec.jl>>[init]
module BoxesSpec
    using Test
    using CarboKitten.Components.Common
    using CarboKitten.Components.Boxes

    @testset "Components/Boxes" begin
        let box = Box{Periodic{2}}(grid_size=(10, 10), phys_scale=2.0u"m"),
            input = Boxes.Input(box=box)
            @test input.box.grid_size == (10, 10)
            @test input.box.phys_scale == 2.0u"m"
        end
    end
end
# ~/~ end