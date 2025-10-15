# ~/~ begin <<docs/src/visualization.md#examples/visualization/production_curve.jl>>[init]

using CairoMakie
using CarboKitten.Visualization: production_curve

save("docs/src/_fig/production_curve.svg", production_curve("data/output/alcap-example.h5"))
# ~/~ end
