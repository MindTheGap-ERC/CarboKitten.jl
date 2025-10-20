
using CarboKitten
using CarboKitten.Components.TimeIntegration: time_axis
using CarboKitten.Export: read_slice, data_export, CSV
using Unitful
using Interpolations
using DelimitedFiles

path = "data/init_topo/example_init_topo.csv"
function initial_topography(path)
    data_matrix, _ = readdlm(path, ',', header=true)

    parsed_matrix = [parse(Float64, strip(replace(s, " m" => ""))) for s in data_matrix]

    return parsed_matrix .* u"m"
	
end

const time = TimeProperties(
		Δt = 200u"yr",
		steps = 12900
	)

function get_SL()

	input_sl = readdlm("data/init_topo/miller_sl.txt", '\t', header=false) 
	time_vector = collect(time_axis(time)) / u"yr" .|> NoUnits
	interpolator = linear_interpolation(time_vector, vec(input_sl))
    return t -> interpolator(ustrip(u"yr", t)) * u"m"
end

function main()

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
		tag = ARGS[1],
		time = time,
		box = Box{Coast}(grid_size = (100, 70), phys_scale = 170.0u"m"),
		facies = facies,
		sea_level = get_SL(),
		initial_topography = initial_topography(),
# Save sections at 4, 8 and 12 km away from the shore
		output=Dict(
            :profile => OutputSpec(slice = (:, 35), write_interval = 1)),
		subsidence_rate = 50.0u"m/Myr",
		insolation = 400.0u"W/m^2",
		
		sediment_buffer_size = 50,
		transport_solver=Val{:forward_euler},
        intertidal_zone=0.0u"m",
		depositional_resolution = 0.5u"m",
		#cementation_time=50.0u"yr",
		disintegration_rate = 100.0u"m/Myr")

	run_model(Model{ALCAP}, input, "data/init_topo/miller_sl.h5")
	header, profile = read_slice("data/init_topo/miller_sl.h5", :profile)
    columns = [profile[i] for i in [23, 47, 71]]
    data_export(
        CSV(
            :stratigraphic_column => "data/init_topo/miller_sl_sc.csv",
            :metadata => "data/init_topo/miller_sl.toml"),
         header,
         columns)
end

main()