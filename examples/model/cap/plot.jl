# ~/~ begin <<docs/src/ca-with-production.md#examples/model/cap/plot.jl>>[init]
#| creates: docs/src/_fig/cap1-summary.png
#| requires: data/output/cap1.h5
#| collect: figures
using GLMakie
using CarboKitten.Visualization

save("docs/src/_fig/cap1-summary.png", summary_plot("data/output/cap1.h5"))
# ~/~ end