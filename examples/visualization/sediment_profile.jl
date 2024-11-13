# ~/~ begin <<docs/src/visualization.md#examples/visualization/sediment_profile.jl>>[init]
#| creates: docs/src/_fig/sediment_profile.png
#| requires: data/output/alcap-example.h5
#| collect: figures

using CairoMakie
using CarboKitten.Export: read_slice
using CarboKitten.Visualization: sediment_profile

save("docs/src/_fig/sediment_profile.png",
    sediment_profile(read_slice("data/output/alcap-example.h5", :, 25)...))
# ~/~ end
