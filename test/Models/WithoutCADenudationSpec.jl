module withoutCADenudationSpec

using Test
using Unitful

using CarboKitten
using CarboKitten.Models: WithoutCADenudation as wD
using CarboKitten.Denudation: EmpiricalDenudation

const FACIES = [
	    wD.Facies(
        production=BenthicProduction(
            maximum_growth_rate = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        initial_sediment = 5.0u"m",
        reactive_surface=10u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.001u"m/yr"
        ),
	    wD.Facies(
        production=BenthicProduction(
            maximum_growth_rate = 400u"m/Myr",
            extinction_coefficient = 0.1u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        initial_sediment = 5.0u"m",
        reactive_surface=10u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.001u"m/yr"
        ),
	    wD.Facies(
        production=BenthicProduction(
            maximum_growth_rate = 100u"m/Myr",
            extinction_coefficient = 0.005u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        initial_sediment = 5.0u"m",
        reactive_surface=10u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.001u"m/yr"
        )]

    const INPUT = wD.Input(
		tag = "withoutCADenudation_test",
		box = Box{Coast}(grid_size=(1, 1), phys_scale=150.0u"m"),
		time = TimeProperties(
			Δt = 200.0u"yr",
			steps = 10),
        output = Dict(
            :full => OutputSpec(write_interval = 1)),
		initial_topography = (x, y) -> 0.0u"m",
		insolation = 400.0u"W/m^2",
		facies = FACIES,
        denudation = EmpiricalDenudation(precip = 800.0u"mm/yr"))

    const OUT = run_model(Model{wD}, INPUT, MemoryOutput(INPUT))
    
    @testset "Models/withoutCADenudation" begin
        @test all(size(OUT.data_volumes[:full].bathymetry) .== (1,1,11))
    end


end