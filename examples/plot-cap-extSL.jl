module Script
using CarboKitten.Visualization
using GLMakie
import OpenScienceFramework as OSF


function main()
    f = Figure()
    plot_crosssection(f[1, 1], "data/ca-extSL.h5")
    save("you-path", f)

end
end

Script.main()
# ~/~ end

 