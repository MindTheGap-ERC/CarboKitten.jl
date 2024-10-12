# ~/~ begin <<docs/src/model-alcap.md#examples/model/alcap/plot.jl>>[init]
#| creates: docs/src/_fig/alcaps-alternative.png
#| requires: data/output/alcap2.h5
#| collect: figures

using CairoMakie
using CarboKitten.Visualization

function main()
  sediment_profile("data/output/alcaps2.h5", 25)
  save("docs/src/_fig/alcaps-alternative.png", fig)
end

main()
# ~/~ end