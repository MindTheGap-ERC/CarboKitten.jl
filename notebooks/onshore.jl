### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ a2023b74-a36c-11ef-0520-b3ac38491c5a
using Pkg; Pkg.activate("../workenv")

# ╔═╡ 02bf7fda-d6d4-4e13-ac30-e4f0c0f5821b
using Revise

# ╔═╡ b1348b07-5778-45e4-96be-be4f401b9431
using CarboKitten

# ╔═╡ 394b7989-ebd4-4aa2-8259-f52c701a2102
using CarboKitten.Models: ALCAP, OnshoreTransport as OT

# ╔═╡ 949e7a15-6c26-4cb5-85cb-249fdd24595a
using GLMakie

# ╔═╡ 26014b0c-b1a6-4c14-b996-5cc98af32a49
using CarboKitten.Visualization: summary_plot

# ╔═╡ 8ce964e6-712c-495d-9379-ed07ca88ad09
using Unitful

# ╔═╡ fdd339db-988a-4b1e-9253-237c38e758cb
using GeometryBasics: Vec2d

# ╔═╡ 5eb0ace5-c35f-4ffa-9081-3fd712960339
using Unitful: μm, yr, s, m, π, cosh, sinh, @u_str

# ╔═╡ 8285516a-650d-4392-a4a4-342bb64b396c
using HDF5

# ╔═╡ fe8114d0-9776-4810-9f57-1a89e86c66dd
# ╠═╡ disabled = true
#=╠═╡
alcap_output = run_model(Model{ALCAP}, ALCAP.Example.INPUT, "../data/output/alcap.h5")
  ╠═╡ =#

# ╔═╡ d3353fee-3199-4a0c-baa7-bea0cb228b13
#=╠═╡
summary_plot(alcap_output)
  ╠═╡ =#

# ╔═╡ ceedbedc-418f-442a-91d1-ded003031593
const VelocityVector = typeof(Vec2d(0, 0)u"m/Myr")

# ╔═╡ 176cc53a-94a7-4a72-b06a-f9a89edd556e
const ShearVector = typeof(Vec2d(0, 0)u"1/Myr")

# ╔═╡ 1ec99563-b18a-4844-b75f-88fa3354296b
function ov(lambda, T, wave_amp, water_depth; calibration_factor=3.e-7)
	function (z)
		# Constants and calculations
		k = 2π / lambda        # Wave number [1/m]
		σ = 2π / T             # Angular frequency [1/s]
		a = wave_amp           # Wave amplitude [m]
		h = water_depth        # Water depth [m]

		# Stokes drift [m/s], scaled by the calibration factor
		u_s = calibration_factor * (1 / 2) * σ * k * a^2 * cosh(2 * k * (z + h)) /
			  sinh(k * h)^2

		# Derivative of Stokes drift wrt depth [1/s], 
		# scaled by the calibration factor
		du_s_dz = calibration_factor * σ * k^2 * a^2 * sinh(2 * k * (z + h)) /
				  sinh(k * h)^2

		# Return tuple for velocity and its gradient
		(VelocityVector(u_s, 0.0u"m/yr"), ShearVector(du_s_dz, 0.0u"yr^-1"))
	end
end

# ╔═╡ 8f95eb18-37a3-4771-a1b7-7302068cf9ee
let
	fig = Figure()
	ax = Axis(fig[1,1])
	x = LinRange(-0.0004π, 0.0004π, 1000)
	#lines!(ax, x, cosh.(4x))
	#lines!(ax, x, sinh.(x).^2)
	#ines!(ax, x, cosh.(4x) ./ sinh.(x).^2)
	lines!(ax, x, (1)./sinh.(x).^2)
	fig
end

# ╔═╡ 057a7123-3480-45dc-a1ed-eaa1c4a24fe9
box = CarboKitten.Box{Coast}(grid_size=(50, 50), phys_scale=300.0u"m")

# ╔═╡ cd45ed7e-a016-4482-b8b7-fa356bd356a8
let
	# x, _ = box_axes(box)
	z = LinRange(0.0, 50.0, 100)u"m"
	f = ov(40.0u"m", 10.0u"s", 1.0u"m", 10.0u"m")
	v = f.(z) |> stack
	getindex.(v[1,:], 1) |> in_units_of(u"m/yr")
	
	fig = Figure()
	ax = Axis(fig[1, 1]) # , limits=(nothing, (0, 100.0)))
	lines!(ax, z |> in_units_of(u"m"), getindex.(v[1,:], 1) |> in_units_of(u"m/yr"))
	fig
end

# ╔═╡ 44727aa5-081a-4c70-803f-47b740dec404
const g = 9.8u"m/s^2"

# ╔═╡ 83faac1c-c99a-465a-ad17-cef6b6046906
v(A, k, w) = A * tanh(k * w) * exp(-k * w)

# ╔═╡ b30257f7-f919-42e3-a88c-b9da652f716b
let
    fig = Figure()
    ax = Axis(fig[1,1], yreversed=true, xlabel="v [m/s]", ylabel="wd [m]")

    wd = LinRange(0, 100, 10000)u"m"
    lines!(ax, v.(1.0u"m/s" * 3.331, 0.02u"1/m", wd) / u"m/s" .|> NoUnits, wd / u"m" .|> NoUnits)

    fig
end

# ╔═╡ ddea7129-ce05-40bd-a9b6-7f1dd184863b
LinRange()

# ╔═╡ e1d2171e-eef4-423a-acf9-1cbb4de2b5bd
ot_output = let	
	function ov(lambda, T, wave_amp, water_depth, cutoff_depth; calibration_factor=3.e-7)
	    function (z)
	        if z > cutoff_depth
	            # Return zero velocity and gradient if depth is greater than cutoff
	            return [(0.0u"m/yr", 0.0u"m/yr"), (0.0u"yr^-1", 0.0u"yr^-1")]
	        end
	
	        # Constants and calculations
	        k = 2π / lambda        # Wave number [1/m]
	        σ = 2π / T             # Angular frequency [1/s]
	        a = wave_amp           # Wave amplitude [m]
	        h = water_depth        # Water depth [m]
	
	        # Stokes drift [m/s], scaled by the calibration factor
	        u_s = calibration_factor * (1 / 2) * σ * k * a^2 * cosh(2 * k * (z + h)) /
	              sinh(k * h)^2
	
	        # Derivative of Stokes drift wrt depth [1/s], 
			# scaled by the calibration factor
	        du_s_dz = calibration_factor * σ * k^2 * a^2 * sinh(2 * k * (z + h)) /
	                  sinh(k * h)^2
	
	        # Return tuple for velocity and its gradient
	        [(u_s, 0.0u"m/yr"), (du_s_dz, 0.0u"yr^-1")]
	    end
	end

	
	FACIES = [
	    OT.Facies(
        maximum_growth_rate = 500u"m/Myr",
        extinction_coefficient = 0.8u"m^-1",
        saturation_intensity = 60u"W/m^2",
		diffusion_coefficient = 10.0u"m/yr",
		onshore_velocity = ov(40.0u"m", 10.0u"s", 1.0u"m", 10.0u"m", 5.0u"m")), 

	    OT.Facies(
        maximum_growth_rate = 400u"m/Myr",
        extinction_coefficient = 0.1u"m^-1",
        saturation_intensity = 60u"W/m^2",
		diffusion_coefficient = 10.0u"m/yr",
	    onshore_velocity = ov(40.0u"m", 10.0u"s", 1.0u"m", 10.0u"m", 5.0u"m")), 

	    OT.Facies(
        maximum_growth_rate = 100u"m/Myr",
        extinction_coefficient = 0.005u"m^-1",
        saturation_intensity = 60u"W/m^2",
		diffusion_coefficient = 10.0u"m/yr",
		onshore_velocity = ov(40.0u"m", 10.0u"s", 1.0u"m", 10.0u"m", 5.0u"m"))
	]

	function sea_level(t)
		10.0u"m" * sin(2π * t / 0.2u"Myr") + 3.0u"m" * sin(2π * t / 0.03u"Myr")
	end
	
	INPUT = OT.Input(
		tag = "ot1",
		box = CarboKitten.Box{Coast}(grid_size=(50, 50), phys_scale=300.0u"m"),
		time = TimeProperties(
			Δt = 100.0u"yr",
			steps = 5000,
			write_interval = 10),
		sea_level = sea_level,
		initial_topography = (x, y) -> - x / 300.0,
		subsidence_rate = 50.0u"m/Myr",
		insolation = 400.0u"W/m^2",
		facies = FACIES,
		depositional_resolution = 0.5u"m",
		sediment_buffer_size = 50,
		disintegration_rate = 50.0u"m/Myr")

	run_model(Model{OT}, INPUT, "ot.h5")
end

# ╔═╡ 11fa7c68-40d0-49d8-b448-9db9626cc3dd
summary_plot(ot_output)

# ╔═╡ e1f24dbb-cd4c-4d0b-be59-c419a0d09297
make_rgb(r, g, b) = let s = r+g+b
	RGBf(r/s, g/s, b/s)
end

# ╔═╡ 29652971-bead-4435-bd2e-9401939c3ae2
h5open(ot_output) do fid
	dep = fid["deposition"][:,:,:,end]
	# fac = getindex.(argmax(dep; dims=1)[1, :, :, :], 1)[:,:,1]
	img = make_rgb.(eachslice(dep, dims=1)...)
	image(img, interpolate=false)
	# heatmap(fac, colormap=cgrad(Makie.wong_colors()[1:3]))
end

# ╔═╡ df49944e-d607-4f4d-8f80-de0cae5a8184


# ╔═╡ Cell order:
# ╠═a2023b74-a36c-11ef-0520-b3ac38491c5a
# ╠═02bf7fda-d6d4-4e13-ac30-e4f0c0f5821b
# ╠═b1348b07-5778-45e4-96be-be4f401b9431
# ╠═394b7989-ebd4-4aa2-8259-f52c701a2102
# ╠═949e7a15-6c26-4cb5-85cb-249fdd24595a
# ╠═26014b0c-b1a6-4c14-b996-5cc98af32a49
# ╠═8ce964e6-712c-495d-9379-ed07ca88ad09
# ╠═fdd339db-988a-4b1e-9253-237c38e758cb
# ╠═fe8114d0-9776-4810-9f57-1a89e86c66dd
# ╠═d3353fee-3199-4a0c-baa7-bea0cb228b13
# ╠═5eb0ace5-c35f-4ffa-9081-3fd712960339
# ╠═ceedbedc-418f-442a-91d1-ded003031593
# ╠═176cc53a-94a7-4a72-b06a-f9a89edd556e
# ╠═1ec99563-b18a-4844-b75f-88fa3354296b
# ╠═8f95eb18-37a3-4771-a1b7-7302068cf9ee
# ╠═057a7123-3480-45dc-a1ed-eaa1c4a24fe9
# ╠═cd45ed7e-a016-4482-b8b7-fa356bd356a8
# ╠═44727aa5-081a-4c70-803f-47b740dec404
# ╠═83faac1c-c99a-465a-ad17-cef6b6046906
# ╠═b30257f7-f919-42e3-a88c-b9da652f716b
# ╠═ddea7129-ce05-40bd-a9b6-7f1dd184863b
# ╠═e1d2171e-eef4-423a-acf9-1cbb4de2b5bd
# ╠═11fa7c68-40d0-49d8-b448-9db9626cc3dd
# ╠═8285516a-650d-4392-a4a4-342bb64b396c
# ╠═29652971-bead-4435-bd2e-9401939c3ae2
# ╠═e1f24dbb-cd4c-4d0b-be59-c419a0d09297
# ╠═df49944e-d607-4f4d-8f80-de0cae5a8184
