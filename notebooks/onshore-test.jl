### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# ╔═╡ 1bd14d3c-f480-11ef-1a4e-e170792388dd
using Pkg; Pkg.activate("../workenv")

# ╔═╡ b5828d77-90b1-476c-976a-52cf5d471515
using Revise

# ╔═╡ 59dfd9ed-b91e-480d-8585-f358274f2b92
using CarboKitten

# ╔═╡ 32d22ca3-a1f5-4796-9b23-df8702374ef0
using CarboKitten.Models: ALCAP

# ╔═╡ 737dee0f-6a96-49e0-80cb-19e727fe34c7
using CarboKitten.Visualization: summary_plot

# ╔═╡ a52fdcac-c5a4-47b2-b198-e82bb4e2b20c
using GLMakie

# ╔═╡ 7afe9f13-c15e-4571-964c-5ffce8a12f00
using Unitful

# ╔═╡ 8b71a88e-5e92-47f2-988a-3ee74ef53631
using Printf

# ╔═╡ e1df1edc-f84f-477e-8631-4856b028d48f
box_axes

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

# ╔═╡ 27281bf9-e806-46d4-967d-6b29ae6f38a6
function plot_erosion2(input, every=40)
	y_idx = 1

	(x, y) = box_axes(input.box)

	fig = Figure()
	ax2 = Axis(fig[1, 1], xlabel="x (km)", ylabel="η (m)")

	state = ALCAP.initial_state(input)

	plot_state() = begin
		t = state.step * input.time.Δt
		η = input.initial_topography.(x, y') .+ state.sediment_height .- input.subsidence_rate * t
		lines!(ax2, x |> in_units_of(u"km"), η[:, y_idx] |> in_units_of(u"m"), label=@sprintf("%.3f Myr", ustrip(t)))
	end

	plot_state()
	run_model(Model{ALCAP}, input, state) do i, _
		if mod(i, every) == 0
			plot_state()
		end
	end

	Legend(fig[1, 2], ax2)

	fig
end

# ╔═╡ a12a9c01-eaa7-4870-a221-6164a5717f94
transport_test_input(;
	initial_topography = (x, y) -> 0.0u"m",
	initial_sediment = (x, y) -> 0.0u"m",
	disintegration_rate = 50.0u"m/Myr",
	subsidence_rate = 50.0u"m/Myr",
	diffusion_coefficient = 0.0u"m/yr",
	wave_velocity = _ -> (Vec2(0.0, 0.0)u"m/yr", Vec2(0.0, 0.0)u"1/yr")) =
	
	ALCAP.Input(
		box = CarboKitten.Box{Coast}(grid_size=(100, 1), phys_scale=150.0u"m"),
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
		depositional_resolution = 1000.0u"m")

# ╔═╡ 84155611-61dd-4da0-8848-fbc44abf314b
let
	function initial_sediment(x, y)
	  if x < 5.0u"km"
	    return 30.0u"m"
	  end
	
	  if x > 10.0u"km" && x < 11.0u"km"
	    return 20.0u"m"
	  end
	
	  return 5.0u"m"
	end

	input = transport_test_input(
		initial_topography = (x, y) -> -30.0u"m",
		initial_sediment = initial_sediment,
		diffusion_coefficient = 10.0u"m/yr")

	plot_erosion2(input, 250)
end

# ╔═╡ c8ef0c9b-7600-4363-99a2-2fa8ab3353d4
md"""
## Translation
"""

# ╔═╡ 38060f13-bf99-4f4a-9b73-1e28d77e8656
let

function gaussian_initial_sediment(x, y)
	exp(-(x-10u"km")^2 / (2 * (0.5u"km")^2)) * 30.0u"m"
end

v_prof(v_max, max_depth, w) = 
	let k = sqrt(0.5) / max_depth,
		A = 3.331 * v_max,
		α = tanh(k * w),
		β = exp(-k * w)
		(A * α * β, -A * k * β * (1 - α - α^2))
	end

v_const(v_max) = _ -> (Vec2(v_max, 0.0u"m/yr"), Vec2(0.0u"1/yr", 0.0u"1/yr"))

input = transport_test_input(
	initial_topography = (x, y)  -> -30.0u"m",
	initial_sediment = gaussian_initial_sediment,
	disintegration_rate = 50000.0u"m/Myr",
	wave_velocity = v_const(-5u"km/Myr")
)
	
plot_erosion2(input, 250)
end

# ╔═╡ b121fa61-c014-4ee0-9db0-b8a0a238b9dd
md"""
## With shear
"""

# ╔═╡ 14504d1a-2803-436b-8fba-6a9a509dfa46
let

v_prof(v_max, max_depth, w) = 
	let k = sqrt(0.5) / max_depth,
		A = 3.331 * v_max,
		α = tanh(k * w),
		β = exp(-k * w)
		(A * α * β, -A * k * β * (1 - α - α^2))
	end

v_const(v_max) = _ -> (Vec2(v_max, 0.0u"m/yr"), Vec2(0.0u"1/yr", 0.0u"1/yr"))

input = transport_test_input(
	initial_topography = (x, y) -> -x / 375.0 - 10u"m",
	initial_sediment = 10.0u"m",
	disintegration_rate = 50.0u"m/Myr",
	diffusion_coefficient = 5.0u"m/yr",
	wave_velocity = v_const(-0.1u"m/yr")
)
	
plot_erosion2(input, 250)
end

# ╔═╡ c9ba2fb2-e29a-4148-b267-33ccb41b66ed
let

v_prof(v_max, max_depth, w) = 
	let k = sqrt(0.5) / max_depth,
		A = 3.331 * v_max,
		α = tanh(k * w),
		β = exp(-k * w)
		(A * α * β, -A * k * β * (1 - α - α^2))
	end

v_const(v_max) = _ -> (Vec2(v_max, 0.0u"m/yr"), Vec2(0.0u"1/yr", 0.0u"1/yr"))

v_prof_par(v_max, max_depth) = w -> let (v, s) = v_prof(v_max, max_depth, w)
	(Vec2(v, 0.0u"m/yr"), Vec2(s, 0.0u"1/yr"))
end
	
input = transport_test_input(
	initial_topography = (x, y) -> -x / 375.0 - 10u"m",
	initial_sediment = 10.0u"m",
	disintegration_rate = 50.0u"m/Myr",
	diffusion_coefficient = 5.0u"m/yr",
	wave_velocity = v_prof_par(-0.5u"m/yr", 20u"m")
)
	
plot_erosion2(input, 250)
end

# ╔═╡ Cell order:
# ╠═1bd14d3c-f480-11ef-1a4e-e170792388dd
# ╠═b5828d77-90b1-476c-976a-52cf5d471515
# ╠═59dfd9ed-b91e-480d-8585-f358274f2b92
# ╠═e1df1edc-f84f-477e-8631-4856b028d48f
# ╠═32d22ca3-a1f5-4796-9b23-df8702374ef0
# ╠═737dee0f-6a96-49e0-80cb-19e727fe34c7
# ╠═a52fdcac-c5a4-47b2-b198-e82bb4e2b20c
# ╠═7afe9f13-c15e-4571-964c-5ffce8a12f00
# ╠═9fe4823d-ae06-46b9-9d10-3ffa03b8fe8b
# ╠═5cf16ac9-14fb-4b4e-a74e-fe2c5708bffd
# ╠═825a57a1-28b1-4995-88c3-486a87118787
# ╠═b8e9bf4b-14a3-4049-8bbc-7db47a064bbe
# ╠═5b99875c-97f6-4fed-bb46-18bc88e10a73
# ╠═4ffdeb5b-b0bc-4932-bcb6-a6d50179d812
# ╠═9ea62335-b2d6-4389-ad26-0f88811f50d2
# ╠═27281bf9-e806-46d4-967d-6b29ae6f38a6
# ╠═a12a9c01-eaa7-4870-a221-6164a5717f94
# ╠═84155611-61dd-4da0-8848-fbc44abf314b
# ╟─c8ef0c9b-7600-4363-99a2-2fa8ab3353d4
# ╠═38060f13-bf99-4f4a-9b73-1e28d77e8656
# ╠═8b71a88e-5e92-47f2-988a-3ee74ef53631
# ╠═b121fa61-c014-4ee0-9db0-b8a0a238b9dd
# ╠═14504d1a-2803-436b-8fba-6a9a509dfa46
# ╠═c9ba2fb2-e29a-4148-b267-33ccb41b66ed
