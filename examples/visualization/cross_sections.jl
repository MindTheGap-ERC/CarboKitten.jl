# ~/~ begin <<docs/src/visualization/profiles.md#examples/visualization/cross_sections.jl>>[init]

module Script

using CairoMakie
using CarboKitten.Export: read_volume
using CarboKitten.Visualization: sediment_profile, sediment_proportion

function main()
    header, vol = read_volume("data/output/alcap-example.h5", :topography)

    nx, ny = size(vol.sediment_thickness)[1:2]
    dip    = vol[:, div(ny, 2) + 1]   # dip section at mid-y
    strike = vol[div(nx, 2) + 1, :]   # strike section at mid-x

    save("docs/src/fig/xsec_dip.png",
         sediment_profile(header, dip))

    save("docs/src/fig/xsec_strike.png",
         sediment_profile(header, strike))

    save("docs/src/fig/xsec_proportion_deposited.png",
         sediment_proportion(header, dip, 1; mode=:deposited))

    save("docs/src/fig/xsec_proportion_preserved.png",
         sediment_proportion(header, dip, 1; mode=:preserved))
end

end

Script.main()
# ~/~ end
