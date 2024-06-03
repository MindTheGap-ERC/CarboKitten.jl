# ~/~ begin <<docs/src/boxes.md#test/ConfigSpec.jl>>[init]
@testset "Config" begin
    using CarboKitten.BoundaryTrait
    using CarboKitten.Config: Box
    using CarboKitten.Vectors

    box = Box{Shelf}(
        grid_size = (100, 50),
        phys_scale = 1.0u"km")
    @test box.phys_size == (x=100000.0, y=50000.0)
end
# ~/~ end