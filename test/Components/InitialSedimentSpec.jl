# ~/~ begin <<docs/src/components/initial_sediment.md#test/Components/InitialSedimentSpec.jl>>[init]
@testset "Components/InitialSediment" begin

using Unitful
using CarboKitten.Components.Common
using CarboKitten.Components: InitialSediment as IS

let box = Box{Periodic}(grid_size=(50, 50), phys_scale=100.0u"m"),
    input = record_type(:Input, IS)(
        time = TimeProperties(steps=0, Δt=1.0u"yr"),
        box = box,
        facies = [ record_type(:Facies, IS)(initial_sediment=30.0u"m") ] ),
    # the initial sediment is pushed by `initialize!`
    state = Common.initial_state(record_type(:State, IS), input)

    @test all(isapprox.(state.sediment_thickness, 30.0u"m"))
    @test all(isapprox.(state.sediment_buffer, 1.0))
end

end
# ~/~ end
