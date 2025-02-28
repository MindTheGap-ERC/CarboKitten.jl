### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 1bd14d3c-f480-11ef-1a4e-e170792388dd
using Pkg; Pkg.activate("../workenv")

# ╔═╡ b5828d77-90b1-476c-976a-52cf5d471515
using Revise

# ╔═╡ 59dfd9ed-b91e-480d-8585-f358274f2b92
using CarboKitten

# ╔═╡ 32d22ca3-a1f5-4796-9b23-df8702374ef0
using CarboKitten.Models: OnshoreTransport as OT, ALCAP

# ╔═╡ 737dee0f-6a96-49e0-80cb-19e727fe34c7
using CarboKitten.Visualization: summary_plot

# ╔═╡ a52fdcac-c5a4-47b2-b198-e82bb4e2b20c
using GLMakie

# ╔═╡ 7afe9f13-c15e-4571-964c-5ffce8a12f00
using Unitful

# ╔═╡ 9fe4823d-ae06-46b9-9d10-3ffa03b8fe8b
sea_level(t) = 10.0u"m" * sin(2π * t / 100.0u"kyr")

# ╔═╡ 5cf16ac9-14fb-4b4e-a74e-fe2c5708bffd
initial_topography(x, y) = -((x - 7.5u"km")^2 + (y - 7.5u"km")^2) / 1000.0u"km"

# ╔═╡ 825a57a1-28b1-4995-88c3-486a87118787


# ╔═╡ b8e9bf4b-14a3-4049-8bbc-7db47a064bbe
md"""
Derivative of $\tanh$ is $1 - \tanh^2$, so:

$$v(w) = A \tanh(k w) \exp(- k w)$$
$$\begin{aligned}
v'(w) &= A k [(1 - \tanh^2(k w)) \exp(- k w) - \tanh(k w) \exp(-k w)]\\
      &= A k \exp(-k w) [1 - \tanh^2(k w) - \tanh(k w)]
\end{aligned}$$
"""

# ╔═╡ 5b99875c-97f6-4fed-bb46-18bc88e10a73
v_prof(v_max, max_depth, w) = 
	let k = sqrt(0.5) / max_depth,
		A = 3.331 * v_max,
		α = tanh(k * w),
		β = exp(-k * w)
		(A * α * β, -A * k * β * (1 - α - α^2))
	end

# ╔═╡ 4ffdeb5b-b0bc-4932-bcb6-a6d50179d812
let w = LinRange(0, 100.0, 1000)u"m"
	f = v_prof.(10.0u"m/yr", 20.0u"m", w)

	v = first.(f)
	s = last.(f)
	
	fig = Figure()
	ax1 = Axis(fig[1, 1], title="transport velocity", yreversed=true, xlabel="velocity [m/yr]", ylabel="depth [m]")
	ax2 = Axis(fig[1, 2], title="transport shear", yreversed=true, xlabel="shear [1/yr]", ylabel="depth [m]")
	lines!(ax1, v / u"m/yr", w / u"m")
	lines!(ax2, s / u"1/yr", w / u"m")
	fig
end

# ╔═╡ 9ea62335-b2d6-4389-ad26-0f88811f50d2
velocity = w -> let (v, s) = v_prof(2.0u"m/yr", 10.0u"m", w)
	(Vec2(-v, 0.0u"m/yr"), Vec2(-s, 0.0u"1/yr"))
end

# ╔═╡ c56cf301-4127-464e-93ae-e5466c94d2fa
facies = [
    OT.Facies(
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=50.0u"m/yr",
		onshore_velocity=velocity),
    OT.Facies(
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=25.0u"m/yr",
		onshore_velocity=velocity),
    OT.Facies(
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=12.5u"m/yr",
		onshore_velocity=velocity)
]

# ╔═╡ c85bb265-2937-4c70-886b-06ecf8d195d1
input = OT.Input(
	tag = "ot-test-small",
	box = CarboKitten.Box{Periodic{2}}((100, 100), 150.0u"m"),
	time = TimeProperties(
		Δt = 100.0u"yr",
		steps = 5000,
		write_interval = 1),
	ca_interval = 10,
	facies = facies,
	sea_level = sea_level,
	initial_topography = initial_topography,
	subsidence_rate = 50.0u"m/Myr",
	insolation = 400.0u"W/m^2",
	sediment_buffer_size = 100,
	depositional_resolution = 0.5u"m",
	disintegration_rate = 0.1u"m/kyr"
)

# ╔═╡ 51c63f8c-e2bd-4605-83bd-3fbc367b1a3c
run_model(Model{ALCAP}, input, "ot-test2.h5")

# ╔═╡ 5a2134a0-dc90-4c9c-b83d-4a5171bdb4e7
summary_plot("ot-test2.h5")

# ╔═╡ 8d5725d6-456b-4637-9e35-df2b4fa1f5ed
# ╠═╡ disabled = true
#=╠═╡
summary_plot("ot-test.h5")
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═1bd14d3c-f480-11ef-1a4e-e170792388dd
# ╠═b5828d77-90b1-476c-976a-52cf5d471515
# ╠═59dfd9ed-b91e-480d-8585-f358274f2b92
# ╠═32d22ca3-a1f5-4796-9b23-df8702374ef0
# ╠═737dee0f-6a96-49e0-80cb-19e727fe34c7
# ╠═a52fdcac-c5a4-47b2-b198-e82bb4e2b20c
# ╠═7afe9f13-c15e-4571-964c-5ffce8a12f00
# ╠═9fe4823d-ae06-46b9-9d10-3ffa03b8fe8b
# ╠═5cf16ac9-14fb-4b4e-a74e-fe2c5708bffd
# ╠═825a57a1-28b1-4995-88c3-486a87118787
# ╟─b8e9bf4b-14a3-4049-8bbc-7db47a064bbe
# ╠═5b99875c-97f6-4fed-bb46-18bc88e10a73
# ╠═4ffdeb5b-b0bc-4932-bcb6-a6d50179d812
# ╠═9ea62335-b2d6-4389-ad26-0f88811f50d2
# ╠═c56cf301-4127-464e-93ae-e5466c94d2fa
# ╠═c85bb265-2937-4c70-886b-06ecf8d195d1
# ╠═51c63f8c-e2bd-4605-83bd-3fbc367b1a3c
# ╠═5a2134a0-dc90-4c9c-b83d-4a5171bdb4e7
# ╠═8d5725d6-456b-4637-9e35-df2b4fa1f5ed
