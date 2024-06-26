# ~/~ begin <<docs/src/ca-prod-with-erosion.md#examples/erosion/plot-caps-osc.jl>>[init]
module Script
    using CarboKitten.Visualization
    using GLMakie

    function main()
        f = Figure()
        plot_crosssection(f[1,1], "data/capes-osc.h5")
	    save("docs/src/fig/capesosc-crosssection.png", f)
    end
end

Script.main()
# ~/~ end