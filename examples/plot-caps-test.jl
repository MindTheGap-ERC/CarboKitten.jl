# ~/~ begin <<docs/src/ca-with-production.md#examples/plot-caps-osc.jl>>[init]
module Script
    using CarboKitten.Visualization
    using GLMakie

    function main()
        f = Figure()
        plot_crosssection(f[1,1], "data/caps-test.h5")
	    save("docs/src/fig/b13-capstest-crosssection.png", f)
    end
end

Script.main()
# ~/~ end