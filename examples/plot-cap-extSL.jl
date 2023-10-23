module Script
using CarboKitten.Visualization
using GLMakie
import OpenScienceFramework as OSF


function main()
    f = Figure()
    plot_crosssection(f[1, 1], "data/ca-extSL.h5")
    save("C:/Users/Liu00141/OneDrive - Universiteit Utrecht/Desktop/work_by_week/external plots/b13-crosssection-ext.png", f)

end
end

Script.main()
# ~/~ end

 