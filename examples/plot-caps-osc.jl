# ~/~ begin <<docs/src/ca-with-production.md#examples/plot-caps-osc.jl>>[init]
#| creates: docs/src/fig/b13-capsosc-crosssection.png
#| requires: data/caps-osc.h5 ext/VisualizationExt.jl
#| collect: figures

module Script
    using CairoMakie
    using GeometryBasics
    using CarboKitten.Visualization

    function main()
        f = Figure()
        plot_crosssection(f[1,1], "data/caps-osc.h5")
     save("docs/src/fig/b13-capsosc-crosssection.png", f)
    end
end

Script.main()
# ~/~ end