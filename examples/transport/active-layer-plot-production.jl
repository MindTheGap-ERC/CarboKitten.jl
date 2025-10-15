# ~/~ begin <<docs/src/active-layer-transport.md#examples/transport/active-layer-plot-production.jl>>[init]

include("active-layer.jl")
using Unitful
using CarboKitten.Config: axes
using CarboKitten.Utility: in_units_of
using CairoMakie
using .ActiveLayer: input

function main()
  (x, y) = axes(input.box)
  η = input.initial_topography.(x, y')
  p = input.production.(x, y')

  fig = Figure()
  ax = Axis3(fig[1,1], xlabel="x (km)", ylabel="y (km)", zlabel="η (m)", azimuth=5π/3)
  surface!(ax, x |> in_units_of(u"km"), y |> in_units_of(u"km"), η |> in_units_of(u"m"), color = p |> in_units_of(u"m/Myr"))
  save("docs/src/_fig/active-layer-production-patch.png", fig)
end

main()
# ~/~ end
