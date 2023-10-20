module Script
using CarboKitten.Visualization
using GLMakie

function main()
    f = Figure()
    plot_crosssection(f[1, 1], "data/ca-extSL.h5")
    save("docs/src/fig/b13-crosssection-ext.png", f)
end
end

Script.main()
# ~/~ end