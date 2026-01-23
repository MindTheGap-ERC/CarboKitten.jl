# ~/~ begin <<docs/src/visualization/topography.md#examples/visualization/glamour_view_single.jl>>[init]
module Script

using GLMakie
using CarboKitten
using CarboKitten.Visualization: glamour_view!
using HDF5

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
	box = CarboKitten.Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m"),
	time = TimeProperties(
		Δt = 200.0u"yr",
		steps = 1000),

    output = Dict(
        :topography => OutputSpec(write_interval = 1000),
        :profile => OutputSpec(slice = (:, 25))),

	sea_level = t -> 4.0u"m" * sin(2π * t / 0.2u"Myr"),
	initial_topography = (x, y) -> - x / 300.0,
	subsidence_rate = 50.0u"m/Myr",
	insolation = 400.0u"W/m^2",
	facies = FACIES)

function main()
    fig = Figure()
    ax = Axis3(fig[1,1])
    out = MemoryOutput(INPUT)
    run_model(Model{CAP}, INPUT, out) 
    glamour_view!(ax, out.header, out.data_volumes[:topography])
    save("docs/src/_fig/glamour_view_single.png", fig)
end

end

Script.main()
# ~/~ end
