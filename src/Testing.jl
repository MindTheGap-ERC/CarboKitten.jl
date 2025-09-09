# ~/~ begin <<docs/src/active-layer-transport.md#src/Testing.jl>>[init]
module Testing

using ..CarboKitten
using ..Models.ALCAP
using Unitful
using GeometryBasics

transport_test_input(;
	initial_topography = (x, y) -> 0.0u"m",
	initial_sediment = (x, y) -> 0.0u"m",
	disintegration_rate = 50.0u"m/Myr",
	subsidence_rate = 0.0u"m/Myr",
	diffusion_coefficient = 0.0u"m/yr",
	wave_velocity = _ -> (Vec2(0.0, 0.0)u"m/yr", Vec2(0.0, 0.0)u"1/yr"),
    intertidal_zone = 0.0u"m") =

	ALCAP.Input(
		box = CarboKitten.Box{Coast}(grid_size=(120, 1), phys_scale=125.0u"m"),
		time = TimeProperties(
			Δt = 0.001u"Myr",
			steps = 1000),
		facies = [ALCAP.Facies(
			initial_sediment = initial_sediment,
			diffusion_coefficient = diffusion_coefficient,
			wave_velocity = wave_velocity,
			maximum_growth_rate = 0.0u"m/Myr",
			extinction_coefficient = 0.8u"m^-1",
			saturation_intensity = 60u"W/m^2"
		)],
		disintegration_rate = disintegration_rate,
		initial_topography = initial_topography,
		insolation = 400.0u"W/m^2",
		sediment_buffer_size = 5,
		depositional_resolution = 1000.0u"m",
		transport_solver = Val{:forward_euler},
        subsidence_rate = subsidence_rate,
        intertidal_zone = intertidal_zone)

end
# ~/~ end