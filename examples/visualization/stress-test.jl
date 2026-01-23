# ~/~ begin <<docs/src/visualization/tests.md#examples/visualization/stress-test.jl>>[init]
using GLMakie
using CarboKitten.Visualization

GLMakie.activate!()

fig_erosion = summary_plot("data/output/alcap-extreme-erosion.h5")
save("docs/src/_fig/alcap-extreme-erosion.png", fig_erosion)

fig_oscillation = summary_plot("data/output/alcap-rapid-oscillation.h5")
save("docs/src/_fig/alcap-rapid-oscillation.png", fig_oscillation)
# ~/~ end
