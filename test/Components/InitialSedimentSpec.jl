# ~/~ begin <<docs/src/components/initial_sediment.md#test/Components/InitialSedimentSpec.jl>>[init]
@testset "Components/InitialSediment" begin

using Unitful
using CarboKitten.Components.Common
using CarboKitten.Components: InitialSediment as IS

function initial_state(input::AbstractInput)
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, IS.n_facies(input), input.box.grid_size...)
    state = IS.State(
        step = 0,
        sediment_height = zeros(IS.Sediment, input.box.grid_size...),
        sediment_buffer = sediment_buffer)
    IS.push_initial_sediment!(input, state)
    return state
end

let box = Box{Periodic}(grid_size=(50, 50), phys_scale=100.0u"m"),
    input = IS.Input(
        time = TimeProperties(steps=0, Δt=1.0u"yr"),
        box = box,
        facies = [ IS.Facies(initial_sediment=30u"m") ] ),
    state = initial_state(input)

    @test all(isapprox.(state.sediment_height, 30.0u"m"))
    @test all(isapprox.(state.sediment_buffer, 1.0))
end

end
# ~/~ end
