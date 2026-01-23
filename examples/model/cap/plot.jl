# ~/~ begin <<docs/src/models/ca-with-production.md#examples/model/cap/plot.jl>>[init]
using GLMakie
using CarboKitten.Visualization

GLMakie.activate!()

save("docs/src/_fig/cap1-summary.png", summary_plot("data/output/cap1.h5"))
# ~/~ end
