# ~/~ begin <<docs/src/ca-with-production.md#examples/plot-caps-osc.jl>>[init]
using CarboKitten.Visualization
using GLMakie

function main()
    f = Figure()
    plot_crosssection(f[1,1], "data/caps-osc.h5")
	save("docs/src/fig/b13-capsosc-crosssection.png", f)
end

main()
# ~/~ end