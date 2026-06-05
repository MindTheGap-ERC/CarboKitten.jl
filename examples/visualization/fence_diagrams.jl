# ~/~ begin <<docs/src/visualization/fence-diagrams.md#examples/visualization/fence_diagrams.jl>>[init]

module Script

using Makie
using Unitful
using CarboKitten
using CarboKitten.Export: read_volume, read_slice
using CarboKitten.Visualization: fence_diagram, fence_diagram!

# -- Example 1: categorical colouring from HDF5 file -------------------------

function from_file_categorical()
    fig = fence_diagram(
        "data/output/alcap-example.h5", :topography;
        x_slices            = [10, 30, 50],
        y_slices            = [2.0u"km", 4.0u"km", 6.0u"km"],
        show_unconformities = 10,
        show_coeval_lines   = true,
        show_sealevel       = true,
        color_by            = :facies)
    save("docs/src/fig/fence_diagram_file_cat.png", fig)
end

# -- Example 2: proportional colouring from HDF5 file -------------------------

function from_file_fraction()
    fig = fence_diagram(
        "data/output/alcap-example.h5", :topography;
        x_slices            = [10, 30, 50],
        y_slices            = [2.0u"km", 4.0u"km", 6.0u"km"],
        show_unconformities = 10,
        show_coeval_lines   = true,
        show_sealevel       = true,
        color_by            = :facies_fraction,
        facies              = 2,
        colormap            = :viridis)
    save("docs/src/fig/fence_diagram_file_fraction.png", fig)
end

# -- Example 3: sequence of DataSlice -----------------------------------------
# Build the slice collection explicitly — any iterable works, including
# results from MemoryOutput or slices read from different files.

function from_slice_sequence()
    header, vol = read_volume("data/output/alcap-example.h5", :topography)
    nx, ny = header.grid_size

    # Pick three dip sections and two strike sections
    slices = [
        vol[:, div(ny, 4)],
        vol[:, div(ny, 2)],
        vol[:, 3 * div(ny, 4)],
        vol[div(nx, 3), :],
        vol[2 * div(nx, 3), :],
    ]

    fig = fence_diagram(header, slices;
        color_by = :facies,
        show_unconformities = true,
        show_coeval_lines   = true)
    save("docs/src/fig/fence_diagram_slice_sequence.png", fig)
end

function main()
    from_file_categorical()
    from_file_fraction()
    from_slice_sequence()
end

end

Script.main()
# ~/~ end
