# ~/~ begin <<docs/src/visualization/profiles.md#examples/visualization/profile_proportion.jl>>[init]

module Script
using CairoMakie
using CarboKitten.Export: read_slice
using CarboKitten.Visualization: sediment_proportion

function main()
    header, data = read_slice("data/output/alcap-example.h5", :profile)
    save("docs/src/_fig/profile_proportion.png",
         sediment_proportion(header, data, 1; mode=:preserved))
end
end

Script.main()
# ~/~ end
