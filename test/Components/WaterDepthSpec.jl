# ~/~ begin <<docs/src/components/waterdepth.md#test/Components/WaterDepthSpec.jl>>[init]
using CarboKitten
import CarboKitten.Components.WaterDepth as WD
using CarboKitten.Components.Common: record_type
using CarboKitten.Components.WaterDepth:
    TimeVariableSubsidenceRate, SpaceTimeVariableSubsidenceRate

@testset "Components/WaterDepth" begin
    input = record_type(:Input, WD)(
        box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
        time = TimeProperties(Δt=1.0u"Myr", steps=10),
        sea_level = t -> 2.0u"m",
        initial_topography = (x, y) -> -10.0u"m",
        subsidence_rate = 5.0u"m/Myr"
    )
    state = WD.initial_state(input)


    @test all(state.bathymetry .== WD.initial_topography(input))

    sub! = WD.subsider(input)
    sub!(state)
    @test all(state.bathymetry.==-15.0u"m")

    wd = WD.water_depth(input)
    @test all(wd(state) .== 17.0u"m")

    @test all(WD.cumulative_subsidence(input) .≈ 50.0u"m")
end

@testset "Components/WaterDepth/variable-subsidence" begin
    box = Box{Periodic{2}}(grid_size=(10, 5), phys_scale=1.0u"m")
    time = TimeProperties(t0=1.0u"Myr", Δt=1.0u"Myr", steps=4)
    x, y = box_axes(box)

    function run_subsider(input)
        state = record_type(:State, WD)(step=0, bathymetry=zeros(typeof(1.0u"m"), box.grid_size...))
        sub! = WD.subsider(input)
        for _ in 1:input.time.steps
            sub!(state)
            state.step += 1
        end
        return state.bathymetry
    end

    @testset "time variable" begin
        # rate is 1, 2, 3, 4 m/Myr at t = 1, 2, 3, 4 Myr
<<<<<<< HEAD
        input = record_type(:Input, WD)(
=======
        input = WD.Input(
>>>>>>> main
            box = box, time = time,
            subsidence_rate = TimeVariableSubsidenceRate(t -> t * 1.0u"m/Myr^2"))
        @test all(run_subsider(input) .≈ -10.0u"m")
        @test all(WD.cumulative_subsidence(input) .≈ 10.0u"m")
    end

    @testset "space-time variable" begin
        f(x, y, t) = x * t * 1.0u"1/Myr^2" + y * 1.0u"1/Myr"
<<<<<<< HEAD
        input = record_type(:Input, WD)(
=======
        input = WD.Input(
>>>>>>> main
            box = box, time = time,
            subsidence_rate = SpaceTimeVariableSubsidenceRate(f))
        expected = x .* 10.0u"Myr^-1" .* 1.0u"Myr" .+ y' .* 4.0u"Myr^-1" .* 1.0u"Myr"
        @test run_subsider(input) ≈ .-expected
        @test WD.cumulative_subsidence(input) ≈ expected
    end

    @testset "constant rate is a special case" begin
<<<<<<< HEAD
        constant = record_type(:Input, WD)(box = box, time = time, subsidence_rate = 3.0u"m/Myr")
        timevar = record_type(:Input, WD)(box = box, time = time,
            subsidence_rate = TimeVariableSubsidenceRate(_ -> 3.0u"m/Myr"))
        spacetimevar = record_type(:Input, WD)(box = box, time = time,
=======
        constant = WD.Input(box = box, time = time, subsidence_rate = 3.0u"m/Myr")
        timevar = WD.Input(box = box, time = time,
            subsidence_rate = TimeVariableSubsidenceRate(_ -> 3.0u"m/Myr"))
        spacetimevar = WD.Input(box = box, time = time,
>>>>>>> main
            subsidence_rate = SpaceTimeVariableSubsidenceRate((_, _, _) -> 3.0u"m/Myr"))
        @test run_subsider(timevar) ≈ run_subsider(constant)
        @test run_subsider(spacetimevar) ≈ run_subsider(constant)
    end
end
# ~/~ end
