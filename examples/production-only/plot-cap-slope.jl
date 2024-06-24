# ~/~ begin <<docs/src/ca-with-production.md#examples/production-only/plot-cap-slope.jl>>[init]
#| creates: docs/src/_fig/b13-crosssection.png
#| requires: data/ca-prod-slope.h5 ext/VisualizationExt.jl
#| collect: figures

module Script
using CairoMakie
using CarboKitten.Visualization

function main()
    f = Figure()
    plot_crosssection(f[1, 1], "data/ca-prod-slope.h5")
    save("docs/src/_fig/b13-crosssection.png", f)
end
end

Script.main()
# ~/~ end