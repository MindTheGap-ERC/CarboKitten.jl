# ~/~ begin <<docs/src/bosscher-1992.md#examples/model/bs92/plot.jl>>[init]
#| creates: docs/src/_fig/bs92-summary.png
#| requires: data/output/bs92.h5
#| collect: figures

using GLMakie
using CarboKitten.Visualization

GLMakie.activate!()

save("docs/src/_fig/bs92-summary.png", summary_plot("data/output/bs92.h5"))
# ~/~ end
