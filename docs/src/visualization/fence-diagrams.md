# Fence Diagrams

The fence-diagram visualization generates multiple cross-sections through the
model domain by specifying x and y slice positions. This provides an intuitive
way to inspect 3D stratigraphic architecture and examine lateral and vertical
variations in facies across the platform.

![Fence diagram — categorical](../_fig/fence_diagram_file_cat.png)
![Fence diagram — proportion](../_fig/fence_diagram_file_fraction.png)

### Examples

Example 1 reproduces the fence diagram above from `alcap-example.h5` with
categorical colouring. Each colour corresponds to one facies.

Example 2 reproduces the same diagram but colours fences by the proportion of
a single selected facies, making it easier to visualise its spatial distribution
and relative abundance.

Example 3 demonstrates the **sequence form** of `fence_diagram!`. Instead of
passing a full `DataVolume`, an iterable of pre-selected `DataSlice` objects is
passed directly. This form is useful for combining slices from multiple HDF5
files, plotting all profiles saved to a `MemoryOutput`, or selecting arbitrary
subsets without loading the full volume.

Example 4 shows the **all-slices HDF5 form**: `fence_diagram(filename)` reads
every `DataSlice` group in the file and plots them all without any manual group
listing. A matching `fence_diagram(output::MemoryOutput)` form works the same
way for in-memory results.

``` {.julia .task file=examples/visualization/fence_diagrams.jl}
#| creates: docs/src/_fig/fence_diagram_file_cat.png
#|          docs/src/_fig/fence_diagram_file_fraction.png
#|          docs/src/_fig/fence_diagram_slice_sequence.png
#| requires: data/output/alcap-example.h5
#| collect: figures

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
    save("docs/src/_fig/fence_diagram_file_cat.png", fig)
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
    save("docs/src/_fig/fence_diagram_file_fraction.png", fig)
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
    save("docs/src/_fig/fence_diagram_slice_sequence.png", fig)
end

function main()
    from_file_categorical()
    from_file_fraction()
    from_slice_sequence()
end

end

Script.main()
```

![Fence diagram — slice sequence](../_fig/fence_diagram_slice_sequence.png)

### Implementation

The implementation is split into three layers:

1. **`fence_plot!`** — renders a single `DataSlice` as a 3D mesh panel.
   Re-uses `explode_quad_vertices` from `SedimentProfile` and delegates height
   computation to `surface_heights` (defined in `Output.Abstract`).

2. **`fence_diagram!(ax, header, slices)`** — core sequence method. Iterates
   over any collection of `DataSlice` objects, applying colour functions and
   decorations to each. Bedrock and sea-level surfaces are not added here
   since they require the full domain extent.

3. **`fence_diagram!(ax, header, data::DataVolume)`** — convenience method.
   Converts `x_slices`/`y_slices` arguments into `DataSlice` objects,
   adds bedrock and sea-level context surfaces, then delegates to the core
   sequence method.

#### `surface_heights`

The height array for a slice (shape `(n_pos, n_t+1)`) used to be computed
inside `FenceDiagrams.jl`. It now lives in `Output.Abstract` as
`surface_heights(header, data::Data{F,D})` and works generically for
`DataColumn`, `DataSlice`, and `DataVolume`. It is also available to
`SedimentProfile.jl` and any future visualization that needs the same quantity.

``` {.julia file=ext/FenceDiagrams.jl}
module FenceDiagram

import CarboKitten.Visualization: fence_diagram, fence_diagram!

using CarboKitten.Visualization
using CarboKitten.Utility: in_units_of
using CarboKitten.Export: Header, Data, DataSlice, DataVolume, read_volume
using CarboKitten.Algorithms: skeleton
using CarboKitten.Output.Abstract: stratigraphic_column, water_depth, surface_heights

import ..SedimentProfile: explode_quad_vertices

using Makie
using GeometryBasics
using Unitful

const Rate   = typeof(1.0u"m/Myr")
const Amount = typeof(1.0u"m")
const Length = typeof(1.0u"m")
const Time   = typeof(1.0u"Myr")

<<fence-helpers>>
<<fence-plot>>
<<fence-decorations>>
<<fence-diagram-sequence>>
<<fence-diagram-volume>>
<<fence-diagram-figure>>

end  # module
```
