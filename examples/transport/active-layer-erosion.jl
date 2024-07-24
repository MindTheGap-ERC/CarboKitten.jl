# ~/~ begin <<docs/src/active-layer-transport.md#examples/transport/active-layer-erosion.jl>>[init]
#| requires: examples/transport/active-layer.jl
#| creates: docs/src/_fig/active-layer-erosion.png
#| collect: figures

module ActiveLayerErosion

include("active-layer.jl")

using Unitful
using CarboKitten.BoundaryTrait: Shelf
using CarboKitten.Config: Box, axes
using CarboKitten.Utility: in_units_of
using CairoMakie

# ~/~ begin <<docs/src/active-layer-transport.md#example-active-layer-erosion>>[init]
function initial_sediment(x, y)
  if x < 5.0u"km"
    return 30.0u"m"
  end

  if x > 10.0u"km" && x < 11.0u"km"
    return 20.0u"m"
  end

  return 5.0u"m"
end

const INPUT = ActiveLayer.Input(
	box                   = Box{Shelf}(grid_size=(100, 1), phys_scale=150.0u"m"),
	Δt                    = 0.001u"Myr",
	t_end                 = 1.0u"Myr",

	bedrock_elevation     = (x, y) -> -30.0u"m",
	initial_sediment      = initial_sediment,
	production            = (x, y) -> 0.0u"m/Myr",

	disintegration_rate   = 50.0u"m/Myr",
	subsidence_rate       = 50.0u"m/Myr",
	diffusion_coefficient = 10000.0u"m")
# ~/~ end

function main(input)
    y_idx = 1
    result = Iterators.map(deepcopy,
  	    Iterators.filter(x -> mod(x[1]-1, 400) == 0, enumerate(ActiveLayer.run_model(input)))) |> collect

	(x, y) = axes(input.box)
	# p = input.production.(x, y')

	fig = Figure(size=(800, 600))
	# ax = Axis3(fig[1:2,1], xlabel="x (km)", ylabel="y (km)", zlabel="η (m)", azimuth=5π/3)
	# surface!(ax, x |> in_units_of(u"km"), y |> in_units_of(u"km"), η |> in_units_of(u"m"))

	ax2 = Axis(fig[1,1], xlabel="x (km)", ylabel="η (m)")

	for r in result
		η = input.bedrock_elevation.(x, y') .+ r[2].sediment .- input.subsidence_rate * r[2].time

		lines!(ax2, x |> in_units_of(u"km"), η[:, y_idx] |> in_units_of(u"m"))
	end

	save("docs/src/_fig/active-layer-erosion.png", fig)
end

end

ActiveLayerErosion.main(ActiveLayerErosion.INPUT)
# ~/~ end