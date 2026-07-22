# ~/~ begin <<docs/src/components/waterdepth.md#test/Components/WaterDepthSpec.jl>>[init]
using CarboKitten
import CarboKitten.Components.WaterDepth as WD
using CarboKitten.Components.WaterDepth: MultiplyRate, Halve, cumulative_subsidence

# Original test from main — unchanged
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
    @test all(state.bathymetry .== -15.0u"m")

    wd = WD.water_depth(input)
    @test all(wd(state) .== 17.0u"m")
end

@testset "Components/WaterDepth/matrix_rate" begin
    nx, ny = 5, 3
    rates  = [10.0u"m/Myr" + 10.0u"m/Myr" * (i - 1) for i in 1:nx, _ in 1:ny]
    input  = WD.Input(
        box = Box{Periodic{2}}(grid_size=(nx, ny), phys_scale=1000.0u"m"),
        time = TimeProperties(Δt=1.0u"Myr", steps=5),
        sea_level = t -> 0.0u"m",
        initial_topography = (x, y) -> -20.0u"m",
        subsidence_rate = rates)

    state = WD._initial_state(input)
    WD.subsider(input)(state)
    for i in 1:nx, j in 1:ny
        @test state.bathymetry[i, j] ≈ -20.0u"m" - rates[i, j] * 1.0u"Myr"
    end
end

@testset "Components/WaterDepth/modifier_halve" begin
    input = WD.Input(
        box = Box{Periodic{2}}(grid_size=(6, 2), phys_scale=1.0u"m"),
        time = TimeProperties(Δt=0.5u"Myr", steps=4),
        sea_level = t -> 0.0u"m",
        initial_topography = (x, y) -> -10.0u"m",
        subsidence_rate = 20.0u"m/Myr",
        subsidence_modifiers = [
            Halve(x_range=(0.0u"m", 3.0u"m"),
                  t_range=(0.0u"Myr", 0.5u"Myr"))])

    state = WD._initial_state(input)
    WD.subsider(input)(state)
    @test state.bathymetry[1, 1] ≈ -15.0u"m"   # halved: 10 m/Myr × 0.5 Myr
    @test state.bathymetry[6, 1] ≈ -20.0u"m"   # full:   20 m/Myr × 0.5 Myr
end

@testset "Components/WaterDepth/cumulative_analytic" begin
    base = fill(10.0u"m/Myr", 3, 3)
    x    = (1.0:3.0) * u"m"
    y    = (1.0:3.0) * u"m"
    cum  = cumulative_subsidence(base, [], x, y, 0.0u"Myr")
    @test all(cum(1.0u"Myr") .≈ 10.0u"m")
    @test all(cum(2.0u"Myr") .≈ 20.0u"m")

    # Halve for first 0.5 Myr: 5×0.5 + 10×0.5 = 7.5 m at t = 1 Myr
    mods = [MultiplyRate(0.5; t_range=(0.0u"Myr", 0.5u"Myr"))]
    cum2 = cumulative_subsidence(base, mods, x, y, 0.0u"Myr")
    @test all(cum2(1.0u"Myr") .≈ 7.5u"m")
end
# ~/~ end
