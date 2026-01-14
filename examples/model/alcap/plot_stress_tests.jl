# Plotting script for stress test scenarios

using GLMakie
using CarboKitten.Visualization

GLMakie.activate!()

fig_erosion = summary_plot("data/output/alcap-extreme-erosion.h5")
save("data/output/alcap-extreme-erosion.png", fig_erosion)

fig_oscillation = summary_plot("data/output/alcap-rapid-oscillation.h5")
save("data/output/alcap-rapid-oscillation.png", fig_oscillation)
