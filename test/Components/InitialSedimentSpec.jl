# ~/~ begin <<docs/src/components/initial_sediment.md#test/Components/InitialSedimentSpec.jl>>[init]
@testset "Components/InitialSediment" begin

using Unitful
using CarboKitten.Components.Common
using CarboKitten.Components.WaterDepth: initial_topography
using CarboKitten.Components.FaciesBase: n_facies
using CarboKitten.Components: InitialSediment as IS

function initial_state(input::AbstractInput)
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...)
    state = record_type(:State, IS)(
        step = 0,
        bathymetry = initial_topography(input),
        sediment_thickness = zeros(IS.Sediment, input.box.grid_size...),
        sediment_buffer = sediment_buffer)
    IS.push_initial_sediment!(input, state)
    return state
end

let box = Box{Periodic}(grid_size=(50, 50), phys_scale=100.0u"m"),
    input = record_type(:Input, IS)(
        time = TimeProperties(steps=0, Δt=1.0u"yr"),
        box = box,
        facies = [ record_type(:Facies, IS)(initial_sediment=30.0u"m") ] ),
    state = initial_state(input)

    @test all(isapprox.(state.sediment_thickness, 30.0u"m"))
    @test all(isapprox.(state.sediment_buffer, 1.0))
end

end
# ~/~ end
