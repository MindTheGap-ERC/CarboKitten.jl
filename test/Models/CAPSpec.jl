# ~/~ begin <<docs/src/models/ca-with-production.md#test/Models/CAPSpec.jl>>[init]
module CAPSpec

using Test
using Unitful

using CarboKitten
using CarboKitten.Models: CAP

    const FACIES = [
	    CAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        production=BenthicProduction(
            maximum_growth_rate = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2")
        ),
	    CAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        production=BenthicProduction(
            maximum_growth_rate = 400u"m/Myr",
            extinction_coefficient = 0.1u"m^-1",
            saturation_intensity = 60u"W/m^2")
        ),
	    CAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        production=BenthicProduction(
            maximum_growth_rate = 100u"m/Myr",
            extinction_coefficient = 0.005u"m^-1",
            saturation_intensity = 60u"W/m^2")
        )]

    const INPUT = CAP.Input(
		tag = "cap_test",
		box = Box{Coast}(grid_size=(1, 1), phys_scale=150.0u"m"),
		time = TimeProperties(
			Δt = 200.0u"yr",
			steps = 10),
        output = Dict(
            :full => OutputSpec(write_interval = 1)),
		initial_topography = (x, y) -> - x / 300.0,
		insolation = 400.0u"W/m^2",
		facies = FACIES)

    const OUT = run_model(Model{CAP}, INPUT, MemoryOutput(INPUT))

    @testset "Models/CAP" begin
        @test all(size(OUT.data_volumes[:full].sediment_thickness) .== (1,1,11))
    end

end

# ~/~ end
