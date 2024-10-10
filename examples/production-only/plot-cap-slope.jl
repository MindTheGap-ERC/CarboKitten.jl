# ~/~ begin <<docs/src/ca-with-production.md#examples/production-only/plot-cap-slope.jl>>[init]
#| creates: docs/src/_fig/b13-crosssection.png
#| requires: data/ca-prod-slope.h5 ext/VisualizationExt.jl
#| collect: figures

module Script
using CairoMakie
using CairoMakie.Export: read_slice
using CarboKitten.Visualization

function main()
    header, data = read_slice("data/ca-prod-slope.h5", 25)
    fig = sediment_profile(header, data)
    save("docs/src/_fig/b13-crosssection.png", fig)
end
end

Script.main()
# ~/~ end