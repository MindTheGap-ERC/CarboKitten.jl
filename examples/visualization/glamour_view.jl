# ~/~ begin <<docs/src/visualization/topography.md#examples/visualization/glamour_view.jl>>[init]

module Script

using GLMakie
using CarboKitten.Export: read_volume
using CarboKitten.Visualization: glamour_view!
using HDF5

function main()
    fig = Figure()
    ax = Axis3(fig[1,1])
    header, volume = read_volume("data/output/cap1.h5", :topography)
    glamour_view!(ax, header, volume)
    save("docs/src/_fig/glamour_view.png", fig)
end

end

Script.main()
# ~/~ end
