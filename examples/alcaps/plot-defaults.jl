# ~/~ begin <<docs/src/model-alcap.md#examples/alcaps/plot-defaults.jl>>[init]
#| requires: ext/VisualizationExt.jl data/alcaps_default.h5
#| creates: docs/src/_fig/alcaps_default_profile.png
#| collect: figures

using CairoMakie
using Statistics
using GeometryBasics
using CarboKitten.Visualization

function main()
  fig = Visualization.sediment_profile("data/alcaps_default.h5", 25)
  save("docs/src/_fig/alcaps_default_profile.png", fig)
end

main()
# ~/~ end
