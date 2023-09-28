# ~/~ begin <<docs/src/ca-with-production.md#examples/plot-cap-slope.jl>>[init]
using CarboKitten.Visualization
using GLMakie

function main()
    f = Figure()
    plot_crosssection(f[1,1], "data/ca-prod-slope.h5")
	save("docs/src/fig/b13-crosssection.png", f)
end

main()
# ~/~ end