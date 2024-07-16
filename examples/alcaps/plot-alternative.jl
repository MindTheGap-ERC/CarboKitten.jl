# ~/~ begin <<docs/src/model-alcap.md#examples/alcaps/plot-alternative.jl>>[init]
#| creates: docs/src/_fig/alcaps-alternative.png
#| requires: data/alcaps2.h5
#| collect: figures

using CairoMakie
using Statistics
using GeometryBasics
using CarboKitten.Visualization

function main()
  fig = Visualization.sediment_profile("data/alcaps2.h5", 25)
  save("docs/src/_fig/alcaps-alternative.png", fig)
end

main()
# ~/~ end