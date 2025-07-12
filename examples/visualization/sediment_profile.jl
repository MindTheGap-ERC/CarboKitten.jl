# ~/~ begin <<docs/src/visualization.md#examples/visualization/sediment_profile.jl>>[init]
#| creates: docs/src/_fig/sediment_profile.png
#| requires: data/output/alcap-example.h5
#| collect: figures

module Script
using CairoMakie
using CarboKitten.Export: read_slice
using CarboKitten.Visualization: sediment_profile

function main()
    save("docs/src/_fig/sediment_profile.png",
        sediment_profile(read_slice("data/output/alcap-example.h5", :slice)...))
end
end

Script.main()
# ~/~ end
