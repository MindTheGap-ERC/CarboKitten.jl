# ~/~ begin <<docs/src/components/waterdepth.md#test/Components/WaterDepthSpec.jl>>[init]
using CarboKitten
import CarboKitten.Components.WaterDepth as WD

@testset "Components/WaterDepth" begin
    input = WD.Input(
        box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
        time = TimeProperties(Δt=1.0u"Myr", steps=10),
        sea_level = t -> 2.0u"m",
        initial_topography = (x, y) -> -10.0u"m",
        subsidence_rate = 5.0u"m/Myr"
    )
    state = WD._initial_state(input)


    @test all(state.bathymetry .== WD.initial_topography(input))

    sub! = WD.subsider(input)
    sub!(state)
    @test all(state.bathymetry.==-15.0u"m")

    wd = WD.water_depth(input)
    @test all(wd(state) .== 17.0u"m")
end
# ~/~ end
