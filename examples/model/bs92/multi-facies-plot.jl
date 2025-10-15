# ~/~ begin <<docs/src/bosscher-1992.md#examples/model/bs92/multi-facies-plot.jl>>[init]

using GLMakie
using CarboKitten.Visualization

GLMakie.activate!()

save("docs/src/_fig/bs92-multi-facies.png", summary_plot("data/output/bs92-multi-facies.h5", wheeler_smooth=(3, 5)))
# ~/~ end
