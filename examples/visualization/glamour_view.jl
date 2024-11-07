# ~/~ begin <<docs/src/visualization.md#examples/visualization/glamour_view.jl>>[init]
#| creates: docs/src/_fig/glamour_view.png
#| requires: data/output/cap1.h5
#| collect: figures

module Script

using GLMakie
using CarboKitten.Visualization: glamour_view!
using HDF5

function main()
    fig = Figure()
    ax = Axis3(fig[1,1])
    h5open("data/output/cap1.h5", "r") do fid
        glamour_view!(ax, fid)
    end
    save("docs/src/_fig/glamour_view.png", fig)
end

end

Script.main()
# ~/~ end
