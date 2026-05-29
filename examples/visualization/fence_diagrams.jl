# ~/~ begin <<docs/src/visualization/fence-diagrams.md#examples/visualization/fence_diagrams.jl>>[init]

module Script
 
using WGLMakie
using Unitful
using CarboKitten
using CarboKitten.Export: read_volume
using CarboKitten.Visualization: fence_diagram, fence_diagram!
 
# -----------------------------------------------------------------------------
# Example 1 — straight from an HDF5 file. Categorical colouring. 
# -----------------------------------------------------------------------------
#| creates: docs/src/_fig/fence_diagram_file_cat.png or docs/src/_fig/fence_diagram_file_inplace_cat.png
#| requires: data/output/alcap-example.h5
#| collect: figures

function from_file()
    fig = fence_diagram(
        "data/output/alcap-example.h5", :topography;
        x_slices = [10, 30, 50],     # grid indices
        y_slices = [2.0u"km", 4.0u"km", 6.0u"km"],       # physical positions work too
        show_unconformities = 10,
        show_coeval_lines   = true,
        show_sealevel       = true,
        color_by = :facies)
    save("docs/src/fig/fence_diagram_file_cat.png", fig)
    return fig
end
 
# Same data, but using the in-place / lower-level form so you can compose with
# other axes in your own figure.
function from_file_inplace()
    header, volume = read_volume("data/output/alcap-example.h5", :topography)
 
    fig = Figure(size = (1400, 900))
    ax  = Axis3(fig[1, 1])
    ax.azimuth   = -π/3
    ax.elevation = π/8
 
    fence_diagram!(ax, header, volume;
        x_slices = [10, 30, 50],
        y_slices = [5.0u"km"],
        show_unconformities = 10,
        show_coeval_lines   = true,
        show_bedrock        = true,
        color_by = :facies)
    save("docs/src/fig/fence_diagram_file_inplace_cat.png", fig)
    return fig
end

# -----------------------------------------------------------------------------
# Example 2 — straight from an HDF5 file. Continuous proportional colouring. 
# -----------------------------------------------------------------------------
 #| creates: docs/src/_fig/fence_diagram_file_fraction.png or docs/src/_fig/fence_diagram_file_inplace_fraction.png
#| requires: data/output/alcap-example.h5
#| collect: figures

function from_file()
    fig = fence_diagram(
        "data/output/alcap-example.h5", :topography;
        x_slices = [10, 30, 50],     # grid indices
        y_slices = [2.0u"km", 4.0u"km", 6.0u"km"],       # physical positions work too
        show_unconformities = 10,
        show_coeval_lines   = true,
        show_sealevel       = true,
        color_by = :facies_fraction,
        facies = 2,
        colormap = :viridis)
    save("docs/src/fig/fence_diagram_file_fraction.png", fig)
    return fig
end
 
# Same data, but using the in-place / lower-level form so you can compose with
# other axes in your own figure.
function from_file_inplace()
    header, volume = read_volume("data/output/alcap-example.h5", :topography)
 
    fig = Figure(size = (1400, 900))
    ax  = Axis3(fig[1, 1])
    ax.azimuth   = -π/3
    ax.elevation = π/8
 
    fence_diagram!(ax, header, volume;
        x_slices = [10, 30, 50],
        y_slices = [5.0u"km"],
        show_unconformities = 10,
        show_coeval_lines   = true,
        show_bedrock        = true,
        color_by = :facies_fraction,
        facies = 2,
        colormap = :viridis)
    save("docs/src/fig/fence_diagram_file_inplace_fraction.png", fig)
    return fig
end
# -----------------------------------------------------------------------------
# Example 2 — from MemoryOutput.
# -----------------------------------------------------------------------------
#| collect: figures

function from_memory(result)
    # `result` : object returned by `run_model(..., MemoryOutput(input))`.
    # Pick whichever volume output was registered in `input.output`.
    header = result.header
    volume = result.data_volumes[:topography]
 
    fig = fence_diagram(header, volume;
        x_slices = [div(header.grid_size[1], 4),
                    div(header.grid_size[1], 2),
                    3 * div(header.grid_size[1], 4)],
        y_slices = [div(header.grid_size[2], 2)],
        show_unconformities = true,
        show_coeval_lines   = true,
        show_bedrock        = true,
        color_by = :facies)
    return fig
end
 
end  # module Script
 
Script.from_file()
# ~/~ end
