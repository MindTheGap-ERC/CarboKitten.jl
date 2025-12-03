# ~/~ begin <<docs/src/initial-topography.md#examples/initial_topography/plot.jl>>[init]
module Plot

using GLMakie
using CarboKitten.Visualization: summary_plot

function plot()
    fig = summary_plot("data/output/initial-topography/mainrun.h5")
    save("docs/src/_fig/initial_topography_example.png", fig)
end

end

Plot.plot()
# ~/~ end
