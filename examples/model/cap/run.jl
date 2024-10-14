# ~/~ begin <<docs/src/ca-with-production.md#examples/model/cap/run.jl>>[init]
#| creates: data/output/cap1.h5
#| requires: src/Model/CAP.jl

module Script

using CarboKitten
using CarboKitten.Model.CAP
using CarboKitten.Components.Common
using Unitful

const PERIOD = 200.0u"kyr"
const AMPLITUDE = 4.0u"m"

const FACIES = [
	    CAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate = 500u"m/Myr",
        extinction_coefficient = 0.8u"m^-1",
        saturation_intensity = 60u"W/m^2"),

	    CAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate = 400u"m/Myr",
        extinction_coefficient = 0.1u"m^-1",
        saturation_intensity = 60u"W/m^2"),

	    CAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate = 100u"m/Myr",
        extinction_coefficient = 0.005u"m^-1",
        saturation_intensity = 60u"W/m^2")
	]

	const INPUT = CAP.Input(
		tag = "cap1",
		box = Common.Box{Shelf}(grid_size=(100, 50), phys_scale=150.0u"m"),
		time = TimeProperties(
			Δt = 200.0u"yr",
			steps = 5000,
			write_interval = 10),
		sea_level = t -> 4.0u"m" * sin(2π * t / 0.2u"Myr"),
		bedrock_elevation = (x, y) -> - x / 300.0,
		subsidence_rate = 50.0u"m/Myr",
		insolation = 400.0u"W/m^2",
		facies = FACIES)

	main() = CarboKitten.run(Model{CAP}, INPUT, "data/output/cap1.h5")
end

Script.main()
# ~/~ end
