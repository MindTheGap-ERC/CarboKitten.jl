### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 492d9a64-32d0-11ef-3f66-05d7171abb64
using Pkg; Pkg.activate("../workenv")

# ╔═╡ ae6f74a0-462b-444d-8851-da5a6cdcc167
using Revise

# ╔═╡ 29523801-3bdc-4288-a3a2-40c4c804ceec
using CarboKitten.Stencil: convolution, stencil

# ╔═╡ 0faa746e-1a91-4a58-91d3-59324a6f5c85
using CarboKitten.Config: Box

# ╔═╡ eec47cec-16b7-412a-ba82-04c2fcd18d4e
using CarboKitten.BoundaryTrait: Shelf

# ╔═╡ 9e29a143-7771-4e45-ba10-57ffe58bab97
using CarboKitten.Utility: in_units_of

# ╔═╡ 8fe58aef-cabf-40c2-a38d-1f2f38d5ee89
using Unitful

# ╔═╡ f46428cc-586b-453a-a66e-e991809d723e
using CairoMakie

# ╔═╡ b5743ea2-4905-4009-950e-1d5784361328
md"""
# Active Layer Transport

The following is inspired on well-known **active layer** approaches in river bed sediment transport. All quantities with subscript $f$ are facies dependent. Sediment is measured in meters of deposited material. $P_f$ is the production of sediment per facies in $m/s$. Further unit calculations would be more readable if we consider the unit of sediment as separate, so for instance it doesn't cancel against $m^2$ in the units of sediment flux. In the implementation, $\nu$ has the units of ${\rm m}$ which is totaly weird. TBC

In a model without transport, we could write

$$\sigma + \sum_f {{\partial \eta_f} \over {\partial t}} = \sum_f P_f,$$

where $\sigma$ is the subsidence rate in $m/s$.

We suppose that loose sediment, either fresh production or disintegrated older sediment, is being transported in a layer on top of the sea bed. The flux in this layer is assumed to be directly proportional to the local slope of the sea bed $| \nabla_x \eta_* |$, where $\eta_* = \sum_f \eta_f$, the sum over all facies contributions.

The active layer now contains a concentration $C_f$ particles of different grain size (for each facies $f$). If needed, $C_f = \alpha_f P_f$ where $\alpha_f$ is some facies parameter determining the fraction of production that is available for transport. The sediment flux is given as,

$${\bf q_f} = -\nu_f C_f {\bf \nabla_x} \eta_*.$$

The following is the mass balance:

$$\sigma + \sum_f {{\partial \eta_f} \over {\partial t}} = -\sum_f {\bf \nabla_x} \cdot {\bf q_f} + \sum_f P_f,$$

In our modelling we keep track of individual contributions per facies over time. Note that in other approaches to active layer transport there would be a factor $1/C_f$. Here we have a different interpretation to what the concentration means: the sediment settles down after transport, such that the concentration has no impact on the change in sediment surface elevation.

Combining these equations, and ignoring subsidence for the moment (which is a global effect and can't be expressed on a per-facies basis), we get a component-wise diffusion equation

$${{\partial \eta_f(x)}\over{\partial t}} = {\bf \nabla_x} \cdot \big[ \nu_f \alpha_f\ P_f(x)\ {\bf \nabla_x} \eta_{*}(x) \big] + P_f(x),$$

In our model we need to solve this equation one time-step each iteration. If we solve this using forward methods, we should be reminded of the CFL limit for diffusion equations (depending on the diffusion constants and grid size we shouldn't pick the time steps too large). Alternatively, for these two-dimensional situations, an implicit approach is feasible. Also we should take care that somehow $\nabla(\nu\alpha P \nabla \eta) + P > 0$. The interpretation being that we can't transport more than we produce, even if there is capacity to do so.

To solve this equation, it is nicer to expand the transport-diffusion term using the product rule, in short notation:

$$\partial_t \eta_f = \nu' \nabla P_f(x) \cdot \nabla \eta(x) + \nu' P_f(x) \nabla^2 \eta(x) + P_f(x),$$

where $\nu' = \nu_f \alpha_f$

So we have a advection component with velocity $\nu' \nabla P_f$ and a diffusion component with a coefficient $\nu' P_f$.
"""

# ╔═╡ 848b9520-5a9e-4a91-9e73-0d11185af7fa
function axes(box::Box)
	y_axis = (0:(box.grid_size[2] - 1)) .* box.phys_scale
	x_axis = (0:(box.grid_size[1] - 1)) .* box.phys_scale
	return x_axis, y_axis
end

# ╔═╡ 3481fdac-7da2-49a0-a439-ef860f430298
md"""
## Test 1

Suppose we have an incline in one direction, as per usual on a coastal shelf. Production is happening in a circular patch in our box, with constant rate. In addition, we'll release the top 1m of sediment for further transport.
"""

# ╔═╡ 21b3b7d6-84b0-40c5-a6da-772314cb0212
@kwdef struct Input
	box
	Δt::typeof(1.0u"Myr")
	t_end::typeof(1.0u"Myr")
	bedrock_elevation
	production
	disintegration_rate::typeof(1.0u"m/Myr")
	subsidence_rate::typeof(1.0u"m/Myr")
	diffusion_coefficient::typeof(1.0u"m")
end

# ╔═╡ 8afb521a-6472-4481-a955-cd3700fcc021
production_patch(center, radius, rate) = function(x, y)
	(pcx, pcy) = center
	(x - pcx)^2 + (y - pcy)^2 < radius^2 ?
		rate :
		0.0u"m/Myr"
end

# ╔═╡ 02f95740-810c-44ac-bdd0-5beb50405718
md"""
Establish a grid of 100x50, 15km on each side, dropping from 0 to 50m depth.
"""

# ╔═╡ 425181e0-2dce-4991-937c-ac450f4bfe14
const input = Input(
	box=Box{Shelf}(grid_size=(100, 50), phys_scale=150.0u"m"),
	Δt=0.001u"Myr",
	t_end=1.0u"Myr",
	
	bedrock_elevation = (x, y) -> -x / 300.0,
	
	production = production_patch(
		(5000.0u"m", 3750.0u"m"),
		2.0u"km",
		50.0u"m/Myr"),

	disintegration_rate = 50.0u"m/Myr",
	subsidence_rate = 50.0u"m/Myr",

	diffusion_coefficient = 10000.0u"m"
)

# ╔═╡ dc2645d0-806e-4f0e-b079-4e6eb8b880fb
let
	(x, y) = axes(input.box)
	η = input.bedrock_elevation.(x, y')
	p = input.production.(x, y')
	
	fig = Figure()
	ax = Axis3(fig[1,1], xlabel="x (km)", ylabel="y (km)", zlabel="η (m)", azimuth=5π/3)
	surface!(ax, x |> in_units_of(u"km"), y |> in_units_of(u"km"), η |> in_units_of(u"m"), color = p |> in_units_of(u"m/Myr"))
	fig
end

# ╔═╡ 158761bd-6e81-49bc-9338-d8c3dfd48704
# Not used after all

begin
function ∂x(box::Box{BT}, ::Type{TIn}) where {BT, TIn}
	unit_in = unit(TIn)
	unit_out = unit_in / unit(box.phys_scale)
	TOut = typeof(0.0*unit_out)
	kernel = [-1/2 0 1/2] ./ (box.phys_scale)
	convolution(TIn, TOut, BT, kernel')
end

function ∂y(box::Box{BT}, ::Type{TIn}) where {BT, TIn}
	unit_in = unit(TIn)
	unit_out = unit_in / unit(box.phys_scale)
	TOut = typeof(0.0*unit_out)
	kernel = [-1/2 0 1/2] ./ (box.phys_scale)
	convolution(TIn, TOut, BT, kernel)
end

function ∂2(box::Box{BT}, ::Type{TIn}) where {BT, TIn}
	unit_in = unit(TIn)
	unit_out = unit_in / unit(box.phys_scale^2)
	TOut = typeof(0.0*unit_out)
	kernel = [0 1 0; 1 -4 1; 0 1 0] ./ (box.phys_scale^2)
	convolution(TIn, TOut, BT, kernel)
end
end

# ╔═╡ d20caab9-2b48-41bd-a7b9-54169d8af821
mutable struct State
	time::typeof(1.0u"Myr")
	sediment::Matrix{typeof(1.0u"m")}
end

# ╔═╡ 9a32b0ec-2650-486c-a437-0dfacedfed17
function initial_state(input)
	State(0.0u"Myr", fill(0.0u"m", input.box.grid_size...))
end

# ╔═╡ 4b944698-7e1b-4513-8c2a-9191049e26ca
struct Frame
	t::typeof(1.0u"Myr")
	δ::Matrix{typeof(1.0u"m/Myr")}
end

# ╔═╡ 5912fc1c-7f9e-45ff-a6e4-e6048b71677d
const DeltaT = typeof(1.0u"m/Myr")

# ╔═╡ 57122d72-631b-48df-8274-97f98fb7eb54
const SedT = typeof(1.0u"m")

# ╔═╡ 51febb6b-aeca-4964-842e-687821845dac
md"""
Just as a reminder:

$$\partial_t \eta_f = \nu' \nabla P_f(x) \cdot \nabla \eta(x) + \nu' P_f(x) \nabla^2 \eta(x) + P_f(x)$$

Below is the kernel encoding a central differencing scheme i.e. `[-1, 0, 1]/(2Δx)` for first derivative and `[0 -1 0; -1 4 -1; 0 -1 0]/Δx^2` for the laplacian.
"""

# ╔═╡ 2b66297c-01da-41d5-9edd-0e46e5b0d03e
function pde_stencil(box::Box{BT}, ν) where {BT}
	# ν = input.diffusion_coefficient
	Δx = box.phys_scale

	function kernel(x)
		adv = ν * ((x[3, 2][1] - x[1, 2][1]) * (x[3, 2][2] - x[1, 2][2]) +
				   (x[2, 3][1] - x[2, 1][1]) * (x[2, 3][2] - x[2, 1][2])) /
				  (2Δx)^2
		
		dif = ν * x[2, 2][2] * (x[3, 2][1] + x[2, 3][1] + x[1, 2][1] + 
								x[2, 1][1] - 4*x[2, 2][1]) / (Δx)^2
		
		prd = x[2, 2][2]

		# return adv + dif + prd
		return max(0.0u"m/Myr", adv + dif + prd)
		# return prd
	end

	stencil(Tuple{SedT, DeltaT}, DeltaT, BT, (3, 3), kernel)
end

# ╔═╡ 4b71e31f-83eb-4611-b2b6-52ee5c9c6c97
function propagator(input)
	δ = Matrix{DeltaT}(undef, input.box.grid_size...)
	x, y = axes(input.box)
	μ0 = input.bedrock_elevation.(x, y')

	function active_layer(state)
		max_erosion = input.disintegration_rate * input.Δt
		erosion = min.(max_erosion, state.sediment)
		state.sediment .-= erosion

		# convert erosion back into a rate
		input.production.(x, y') .+ erosion ./ input.Δt
	end

	stc = pde_stencil(input.box, input.diffusion_coefficient)
	apply_pde(μ::Matrix{SedT}, p::Matrix{DeltaT}) = stc(tuple.(μ, p), δ)
	
	function (state)
		p = active_layer(state)
		apply_pde(state.sediment .+ μ0, p)
		return Frame(state.time, δ)
	end
end

# ╔═╡ 7d23b253-47b3-46a2-8f5f-328e91e23e3e
function run_model(input)
	state = initial_state(input)
	prop = propagator(input)
	
	Channel{State}() do ch
		while state.time < input.t_end
			Δ = prop(state)
			state.sediment .+= Δ.δ .* input.Δt
			state.time += input.Δt
			put!(ch, state)
		end
	end
end

# ╔═╡ 7b955eb3-6448-45ac-90da-a81ef320fae4
md"""
We run the model with 1000 time steps but only inspect one in every 100.
"""

# ╔═╡ c4ce2e66-6fd0-40a9-937c-486941410440
result = Iterators.map(deepcopy,
	Iterators.filter(x -> mod(x[1], 100) == 0, enumerate(run_model(input)))) |> collect

# ╔═╡ 0f16d940-c81b-4d8d-a130-a2f4f8e7c723
let
	(x, y) = axes(input.box)
	η = input.bedrock_elevation.(x, y') .+ result[10][2].sediment .- input.subsidence_rate * result[10][2].time
	# p = input.production.(x, y')
	
	fig = Figure(size=(800, 1000))
	ax = Axis3(fig[1:2,1], xlabel="x (km)", ylabel="y (km)", zlabel="η (m)", azimuth=5π/3)
	surface!(ax, x |> in_units_of(u"km"), y |> in_units_of(u"km"), η |> in_units_of(u"m"))

	ax2 = Axis(fig[3,1], xlabel="x (km)", ylabel="η (m)")

	for i in 1:10
		η = input.bedrock_elevation.(x, y') .+ result[i][2].sediment .- input.subsidence_rate * result[i][2].time

		lines!(ax2, x |> in_units_of(u"km"), η[:, 25] |> in_units_of(u"m"))
	end
	fig
end

# ╔═╡ 3efaf9e5-2a64-43a5-bfb7-5754a715f2a8
md"""
Note in the bottom figure, due to sedimentation not keeping up with subsidence, the lines go down in time. We see the sediment transport being favoured to downslope areas, which is what we want. This effect could be made more extreme by increasing the erosion rate.
"""

# ╔═╡ Cell order:
# ╠═492d9a64-32d0-11ef-3f66-05d7171abb64
# ╠═ae6f74a0-462b-444d-8851-da5a6cdcc167
# ╟─b5743ea2-4905-4009-950e-1d5784361328
# ╠═29523801-3bdc-4288-a3a2-40c4c804ceec
# ╠═0faa746e-1a91-4a58-91d3-59324a6f5c85
# ╠═eec47cec-16b7-412a-ba82-04c2fcd18d4e
# ╠═9e29a143-7771-4e45-ba10-57ffe58bab97
# ╠═8fe58aef-cabf-40c2-a38d-1f2f38d5ee89
# ╠═f46428cc-586b-453a-a66e-e991809d723e
# ╠═848b9520-5a9e-4a91-9e73-0d11185af7fa
# ╟─3481fdac-7da2-49a0-a439-ef860f430298
# ╠═21b3b7d6-84b0-40c5-a6da-772314cb0212
# ╠═8afb521a-6472-4481-a955-cd3700fcc021
# ╟─02f95740-810c-44ac-bdd0-5beb50405718
# ╠═425181e0-2dce-4991-937c-ac450f4bfe14
# ╟─dc2645d0-806e-4f0e-b079-4e6eb8b880fb
# ╟─158761bd-6e81-49bc-9338-d8c3dfd48704
# ╠═d20caab9-2b48-41bd-a7b9-54169d8af821
# ╠═9a32b0ec-2650-486c-a437-0dfacedfed17
# ╠═4b944698-7e1b-4513-8c2a-9191049e26ca
# ╠═5912fc1c-7f9e-45ff-a6e4-e6048b71677d
# ╠═57122d72-631b-48df-8274-97f98fb7eb54
# ╟─51febb6b-aeca-4964-842e-687821845dac
# ╠═2b66297c-01da-41d5-9edd-0e46e5b0d03e
# ╠═4b71e31f-83eb-4611-b2b6-52ee5c9c6c97
# ╠═7d23b253-47b3-46a2-8f5f-328e91e23e3e
# ╟─7b955eb3-6448-45ac-90da-a81ef320fae4
# ╠═c4ce2e66-6fd0-40a9-937c-486941410440
# ╟─0f16d940-c81b-4d8d-a130-a2f4f8e7c723
# ╟─3efaf9e5-2a64-43a5-bfb7-5754a715f2a8
