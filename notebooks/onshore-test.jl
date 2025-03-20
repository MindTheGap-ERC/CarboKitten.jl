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

# ╔═╡ ed0e7604-c18f-413d-8386-26d2ec637168
module TransportTest

using Unitful
using CarboKitten: Box, box_axes, Boundary
using CarboKitten.Components.Common
using CarboKitten.BoundaryTrait: Shelf
using CarboKitten.Utility: in_units_of
using CarboKitten.Stencil: stencil!
using CarboKitten.Transport.Advection: transport
using CarboKitten.Transport.Solvers: runge_kutta_4

const Amount = typeof(1.0u"m")
const Length = typeof(1.0u"m")
const Rate = typeof(1.0u"m/Myr")

@kwdef struct Input
    box
    Δt::typeof(1.0u"Myr")
    t_end::typeof(1.0u"Myr")
    initial_topography   # function (x::u"m", y::u"m") -> u"m"
    initial_sediment    # function (x::u"m", y::u"m") -> u"m"
    production          # function (x::u"m", y::u"m") -> u"m/s"
    disintegration_rate::typeof(1.0u"m/Myr")
    subsidence_rate::typeof(1.0u"m/Myr")
    diffusivity::typeof(1.0u"m/yr")
	wave_transport
	transport_substeps::Int = 1
	transport_solver = nothing
end

mutable struct State
    time::typeof(1.0u"Myr")
    sediment::Matrix{typeof(1.0u"m")}
end

function initial_state(input)
    x, y = box_axes(input.box)
    State(0.0u"Myr", input.initial_sediment.(x, y'))
end

struct Frame
    t::typeof(1.0u"Myr")
    δ::Matrix{Amount}
end

function propagator(input)
    x, y = box_axes(input.box)
    μ0 = input.initial_topography.(x, y')
    box = input.box
    Δt = input.Δt
    disintegration_rate = input.disintegration_rate
    production = input.production
    d = input.diffusivity
	v = input.wave_transport
	σ = input.subsidence_rate
	solver = if input.transport_solver === nothing
		runge_kutta_4(Amount, box)
	else
		input.transport_solver
	end

    function active_layer(state)
        max_amount = disintegration_rate * Δt
        amount = min.(max_amount, state.sediment)
        state.sediment .-= amount

        production.(x, y') * Δt .+ amount
    end

    function (state)
        p = active_layer(state)
		water_depth = (σ * state.time) .- (μ0 .+ state.sediment)

		for j in 1:input.transport_substeps
			solver(
				(C, _) -> transport(box, d, v, p, water_depth),
				p, state.time, Δt / input.transport_substeps)
        end

        return Frame(state.time, p)
    end
end

function run_model(input)
    state = initial_state(input)
    prop = propagator(input)

    Channel{State}() do ch
        while state.time < input.t_end
            Δ = prop(state)
            state.sediment .+= Δ.δ
            state.time += input.Δt
            put!(ch, state)
        end
    end
end

end

# ╔═╡ 38060f13-bf99-4f4a-9b73-1e28d77e8656
module Erosion

using Unitful
using ..TransportTest: Input
using CarboKitten: Box, Coast

function initial_sediment(x, y)
  if x < 5.0u"km"
    return 30.0u"m"
  end

  if x > 10.0u"km" && x < 11.0u"km"
    return 20.0u"m"
  end

  return 5.0u"m"
end

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

v_const(v_max) = _ -> ((v_max, 0.0u"m/yr"), (0.0u"1/yr", 0.0u"1/yr"))

const INPUT = Input(
	box                   = Box{Coast}(grid_size=(100, 1), phys_scale=150.0u"m"),
    Δt                    = 0.01u"Myr",
    t_end                 = 1.0u"Myr",

    initial_topography     = (x, y) -> -30.0u"m",
    initial_sediment      = gaussian_initial_sediment,
    production            = (x, y) -> 0.0u"m/Myr",

    disintegration_rate   = 5000.0u"m/Myr",
    subsidence_rate       = 50.0u"m/Myr",
    diffusivity           = 0.0u"m/yr",
	wave_transport        = v_const(-0.005u"m/yr"),
	transport_substeps    = 10
		# w -> let (v, s) = v_prof(-5u"m/yr", 20.0u"m", w)
		# 	((v, 0.0u"m/yr"), (s, 0.0u"1/yr"))
		# end
)

end

# ╔═╡ cac22975-0fa6-4bba-9d2e-7274f18c9799
using ProfileCanvas

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
# ╠═╡ disabled = true
#=╠═╡
run_model(Model{ALCAP}, input, "ot-test2.h5")
  ╠═╡ =#

# ╔═╡ 5a2134a0-dc90-4c9c-b83d-4a5171bdb4e7
# ╠═╡ disabled = true
#=╠═╡
summary_plot("ot-test2.h5")
  ╠═╡ =#

# ╔═╡ 8d5725d6-456b-4637-9e35-df2b4fa1f5ed
# ╠═╡ disabled = true
#=╠═╡
summary_plot("ot-test.h5")
  ╠═╡ =#

# ╔═╡ d77eb785-61f8-43ad-a066-3242f1089cfe


# ╔═╡ c8ef0c9b-7600-4363-99a2-2fa8ab3353d4
md"""
## Erosion
"""

# ╔═╡ 91c5858e-e975-4bcd-b157-010f884a1eac
150.0u"m"/0.0001u"Myr" |> u"m/yr"

# ╔═╡ bc4c03ef-5d53-4a99-8aca-b5cfef07ecea
function plot_erosion(input)
	y_idx = 1
	result = Iterators.map(deepcopy,
		Iterators.filter(x -> mod(x[1]-1, 40) == 0, 
		enumerate(TransportTest.run_model(input)))) |> collect

	(x, y) = box_axes(input.box)

	fig = Figure()
	ax2 = Axis(fig[1, 1], xlabel="x (km)", ylabel="η (m)")

	for r in result
		η = input.initial_topography.(x, y') .+ r[2].sediment .- input.subsidence_rate * r[2].time
		lines!(ax2, x |> in_units_of(u"km"), η[:, y_idx] |> in_units_of(u"m"))
	end

	fig
end

# ╔═╡ bb6c2022-b668-4733-8799-cbe7c8f98231
plot_erosion(Erosion.INPUT)

# ╔═╡ b121fa61-c014-4ee0-9db0-b8a0a238b9dd


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
# ╟─c56cf301-4127-464e-93ae-e5466c94d2fa
# ╟─c85bb265-2937-4c70-886b-06ecf8d195d1
# ╠═51c63f8c-e2bd-4605-83bd-3fbc367b1a3c
# ╠═5a2134a0-dc90-4c9c-b83d-4a5171bdb4e7
# ╠═8d5725d6-456b-4637-9e35-df2b4fa1f5ed
# ╠═ed0e7604-c18f-413d-8386-26d2ec637168
# ╠═d77eb785-61f8-43ad-a066-3242f1089cfe
# ╟─c8ef0c9b-7600-4363-99a2-2fa8ab3353d4
# ╠═38060f13-bf99-4f4a-9b73-1e28d77e8656
# ╠═91c5858e-e975-4bcd-b157-010f884a1eac
# ╠═bc4c03ef-5d53-4a99-8aca-b5cfef07ecea
# ╠═bb6c2022-b668-4733-8799-cbe7c8f98231
# ╠═cac22975-0fa6-4bba-9d2e-7274f18c9799
# ╠═b121fa61-c014-4ee0-9db0-b8a0a238b9dd
