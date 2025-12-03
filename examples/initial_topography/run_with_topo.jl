module MainRun

using GLMakie
using CarboKitten
using CarboKitten.Visualization: summary_plot
using CSV
using Tables

include("get_init_topo.jl")

const PATH = Prerun.PATH

function load_initial_topography()
    (CSV.File(Prerun.DATAFILE) |> Tables.matrix) * u"m"
end

function run()
	facies = [
		ALCAP.Facies(
			maximum_growth_rate = 500.0u"m/Myr",
			extinction_coefficient = 0.8u"m^-1",
			saturation_intensity = 60.0u"W/m^2",
			diffusion_coefficient = 50.0u"m/yr"
		),
		ALCAP.Facies(
			maximum_growth_rate = 400.0u"m/Myr",
			extinction_coefficient = 0.1u"m^-1",
			saturation_intensity = 60.0u"W/m^2",
			diffusion_coefficient = 25.0u"m/yr"
		),
		ALCAP.Facies(
			maximum_growth_rate = 100.0u"m/Myr",
			extinction_coefficient = 0.005u"m^-1",
			saturation_intensity = 60.0u"W/m^2",
			diffusion_coefficient = 12.5u"m/yr"
		)]

	input = ALCAP.Input(
		tag = "mainrun",
		time = time = TimeProperties(
			Δt = 200u"yr",
			steps = 5000
		),
		box = Prerun.INPUT.box,
		facies = facies,
		sea_level = t -> 0.0u"m",
		initial_topography = load_initial_topography(),
		output=Dict(
			:topography => OutputSpec(write_interval=100),
            :profile => OutputSpec(slice = (:, 35))
		),
		subsidence_rate = 50.0u"m/Myr",
		insolation = 400.0u"W/m^2",
		
		transport_solver = Val{:forward_euler},
		sediment_buffer_size = 50,
		depositional_resolution = 0.5u"m",
		cementation_time = 50.0u"yr",
		disintegration_rate = 100.0u"m/Myr")

	run_model(Model{ALCAP}, input, joinpath(PATH, "mainrun.h5"))
end

function plot(main_out)
	fig = summary_plot(main_out)
	save("docs/src/_fig/initial_topography_example.png", fig)
end

end  # module MainRun
