# ~/~ begin <<docs/src/visualization/map-view.md#examples/visualization/map_view.jl>>[init]
module Script

using WGLMakie
using Unitful
using CarboKitten
using CarboKitten.Export: read_volume
using CarboKitten.Visualization: map_view, map_view!

# -----------------------------------------------------------------------------
# Example 1 —  from an HDF5 file. 
# -----------------------------------------------------------------------------
function from_file()
    fig = map_view(
        "data/output/alcap-example.h5", :topography;
        times = [0.2u"Myr", 0.5u"Myr", 1.0u"Myr"],   # mix of physical times…
        show = :preserved,
        show_shoreline = true,
        layout = :row)
    save("docs/src/_fig/map_view_file.png", fig)
    return fig
end

# Same data, in-place form so you can drop a single map into your own figure.
function from_file_inplace()
    header, volume = read_volume("data/output/alcap-example.h5", :topography)

    fig = Figure(size = (900, 700))
    ax  = Axis(fig[1, 1])
    hm  = map_view!(ax, header, volume;
        time = 0.5u"Myr",
        show = :both,
        show_shoreline = true)

    n_facies = size(volume.production, 1)
    Colorbar(fig[1, 2], hm; ticks = 1:n_facies, label = "dominant facies")
    save("docs/src/_fig/map_view_file_inplace.png", fig)
    return fig
end

# -----------------------------------------------------------------------------
# Example 2 — from MemoryOutput. 
# -----------------------------------------------------------------------------
function from_memory(result)
    # `result` is whatever was returned by `run_model(..., MemoryOutput(input))`.
    # Pick the volume key registered in `input.output`.
    header = result.header
    volume = result.data_volumes[:topography]

    # Three evenly-spaced frames through the run, given as integer indices.
    n_frames = length(header.axes.t[1:volume.write_interval:end])
    idx_pick = [div(n_frames, 4), div(n_frames, 2), n_frames]

    fig = map_view(header, volume;
        times = idx_pick,
        show = :preserved,
        layout = :row)
    return fig
end

end  # module Script

Script.from_file()
# ~/~ end
