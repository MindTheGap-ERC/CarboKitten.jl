# ~/~ begin <<docs/src/ca-with-production.md#examples/production-only/plot-no-ca-slope.jl>>[init]
#| creates: docs/src/fig/no-ca-slope.png
#| requires: data/no-ca-slope.h5 ext/VisualizationExt.jl
#| collect: figures

module Script
    using CairoMakie
    using GeometryBasics
    using CarboKitten.Visualization

    function main()
        f = Figure()
        plot_crosssection(f[1,1], "data/no-ca-slope.h5")
	    save("docs/src/fig/no-ca-slope.png", f)
    end
end

Script.main()
# ~/~ end