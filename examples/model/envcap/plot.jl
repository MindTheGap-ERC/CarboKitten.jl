# ~/~ begin <<docs/src/visualization/envcap_plotting.md#examples/model/envcap/plot.jl>>[init]

include("../../../ext/EnvCAPPlotting.jl")

using .EnvCAPPlotting

EnvCAPPlotting.plot_envcap_run_outputs(
    z_index = 41,
    dz = 0.5,
)
# ~/~ end
