# ~/~ begin <<docs/src/active-layer-transport.md#examples/transport/active-layer.jl>>[init]
module ActiveLayer

using CarboKitten.Stencil: convolution, stencil
using CarboKitten.Config: Box, axes
using CarboKitten.BoundaryTrait: Shelf
using CarboKitten.Utility: in_units_of
using Unitful

# ~/~ begin <<docs/src/active-layer-transport.md#example-active-layer>>[init]
using CarboKitten.Stencil: convolution, stencil
using CarboKitten.Config: Box, axes
using CarboKitten.BoundaryTrait: Shelf
using CarboKitten.Utility: in_units_of
using Unitful
# ~/~ end
# ~/~ begin <<docs/src/active-layer-transport.md#example-active-layer>>[1]
@kwdef struct Input
	box
	Δt::typeof(1.0u"Myr")
	t_end::typeof(1.0u"Myr")
	bedrock_elevation   # function (x::u"m", y::u"m") -> u"m"
	production          # function (x::u"m", y::u"m") -> u"m/s"
	disintegration_rate::typeof(1.0u"m/Myr")
	subsidence_rate::typeof(1.0u"m/Myr")
	diffusion_coefficient::typeof(1.0u"m")
end
# ~/~ end
# ~/~ begin <<docs/src/active-layer-transport.md#example-active-layer>>[2]
production_patch(center, radius, rate) = function(x, y)
	(pcx, pcy) = center
	(x - pcx)^2 + (y - pcy)^2 < radius^2 ?
		rate :
		0.0u"m/Myr"
end

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
# ~/~ end
# ~/~ begin <<docs/src/active-layer-transport.md#example-active-layer>>[3]
const DeltaT = typeof(1.0u"m/Myr")
const SedT = typeof(1.0u"m")

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
# ~/~ end
# ~/~ begin <<docs/src/active-layer-transport.md#example-active-layer>>[4]
mutable struct State
	time::typeof(1.0u"Myr")
	sediment::Matrix{typeof(1.0u"m")}
end

function initial_state(input)
	State(0.0u"Myr", fill(0.0u"m", input.box.grid_size...))
end

struct Frame
	t::typeof(1.0u"Myr")
	δ::Matrix{typeof(1.0u"m/Myr")}
end

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
# ~/~ end

end
# ~/~ end