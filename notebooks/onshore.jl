### A Pluto.jl notebook ###
# v0.20.4

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

# ╔═╡ e1d2171e-eef4-423a-acf9-1cbb4de2b5bd
ot_output = let	
	function ov(λ, T, wave_amp, water_depth; cal_fac=1.0e-7)
	    function (z)		

	        # Constants and calculations
	        k = 2π / λ        # Wave number [1/m]
	        σ = 2π / T             # Angular frequency [1/s]
	        a = wave_amp           # Wave amplitude [m]
	        h = water_depth        # Water depth [m]
	
	        # Stokes drift [m/s], scaled by the calibration factor
	        # u_s = -calibration_factor * (1 / 2) * σ * k * a^2 * cosh(2 * k * (z + h))/                sinh(k * h)^2

			u_s = - cal_fac * 4 * π^2 * a^2 * exp(-4 * π * z / λ) / (λ * T)
	
	        # Derivative of Stokes drift wrt depth [1/s], 
			# # scaled by the calibration factor
	  #       du_s_dz = -calibration_factor * σ * k^2 * a^2 * sinh(2 * k * (z + h)) /
	  #                  sinh(k * h)^2

			du_s_dz=  cal_fac * 4^2 * π^3 * a^2 * exp(-4 * π * z / λ) / (λ^2 * T)
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
		onshore_velocity = ov(40.0u"m", 10.0u"s", 1.0u"m", 10.0u"m")), 

	    OT.Facies(
        maximum_growth_rate = 400u"m/Myr",
        extinction_coefficient = 0.1u"m^-1",
        saturation_intensity = 60u"W/m^2",
		diffusion_coefficient = 10.0u"m/yr",
	    onshore_velocity = ov(40.0u"m", 10.0u"s", 1.0u"m", 10.0u"m")), 

	    OT.Facies(
        maximum_growth_rate = 100u"m/Myr",
        extinction_coefficient = 0.005u"m^-1",
        saturation_intensity = 60u"W/m^2",
		diffusion_coefficient = 10.0u"m/yr",
		onshore_velocity = ov(40.0u"m", 10.0u"s", 1.0u"m", 10.0u"m"))
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
# ╠═fe8114d0-9776-4810-9f57-1a89e86c66dd
# ╠═d3353fee-3199-4a0c-baa7-bea0cb228b13
# ╠═5eb0ace5-c35f-4ffa-9081-3fd712960339
# ╠═e1d2171e-eef4-423a-acf9-1cbb4de2b5bd
# ╠═11fa7c68-40d0-49d8-b448-9db9626cc3dd
# ╠═8285516a-650d-4392-a4a4-342bb64b396c
# ╠═29652971-bead-4435-bd2e-9401939c3ae2
# ╠═e1f24dbb-cd4c-4d0b-be59-c419a0d09297
# ╠═df49944e-d607-4f4d-8f80-de0cae5a8184
