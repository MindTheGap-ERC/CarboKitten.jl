Post-processing and Analysis
============================

After CarboKitten has run your model, we then need to do post-processing, analysis and visualization. In the input specification there is a section that determines the shape of all output written to a file (or memory). For example, the following output specification,

```julia
...
    output = Dict(
          :topography => OutputSpec(write_interval = 200),
          :profile => OutputSpec(slice = (:, 75)),
          ("col$(i)" |> Symbol => OutputSpec(slice = (i * 15, 50)) for i in 1:3)...,
...
```

Will create an output with `:topography`, `:profile`, `:col1`, `:col2` and `:col3` data sets.

Loading data
------------

We can use the `load`, `load_volume`, `load_slice` and `load_column` functions to load these from any output.

The `load` function takes either a filename (i.e. `<: AbstractString`) or an output (i.e. `<: AbstractOutput`) and returns a `Bundle`. The bundle can be further queried for any present data sets, which are not loaded into memory until requested.

The `load_volume`, `load_slice` and `load_column` functions take the combination of a filename, `AbstractOutput` or `Bundle` and a symbol to load header and data into a `Pack`. In the case of loading from a `Bundle`, the `Header` information will be shared between them.

A `Pack` consists of a `Header` and `Data` object.

Post-processing
---------------

The following methods are defined on `Pack`s and/or `Header`/`Data` pairs.

| Method | Description |
|---|---|
| `bathymetry` | Subsidence corrected height of sea floor. |
| `disintegration` | Sediments removed at each time step. |
| `production` | Sediments produced at each time step. |
| `deposition` | Sediments deposited at each time step (after transport). |
| `active_layer` | Instantaneous contents of the active layer at each time step (if saved, returns `nothing` otherwise). |
| `sediment_thickness` | Computes the sediment thickness over time, i.e. cumulative sum of sedimentation. |
| `water_depth` | Historic water depth over time. |
| `stratigraphic_column` | Preserved deposition over time. |
| `surface_heights` | Computes the height of preserved surfaces at $t_{\rm end}$. |
| `net_deposition` | Deposition minus disintegration at each time step. |

Some of these may be cached in the `Data` object for efficiency.

The following methods give information on the layout of a `Pack`, without needing to access the members of `Header` or `Data`. Methods marked with (H) are also defined on a `Header`, in which case they refer to the entire model domain and every time step.

| Method | Description |
|---|---|
| `header` | The `Header` of a `Pack`. |
| `n_facies` (H) | Number of facies. |
| `slice` | The slice specification, e.g. `(:, 25)`. |
| `write_interval` | Number of model time steps between writes. |
| `write_time_step` | Time between two writes. |
| `time_axis` (H) | Times at which data was written. |
| `x_axis`, `y_axis` (H) | Coordinates covered by the data; a single coordinate along a sliced axis. |
| `section_axis` | For a slice, the coordinates along which the section runs. |
| `sea_level` (H) | Sea level at each write. |
| `total_subsidence` (H only) | Total subsidence over the run. |
| `bedrock` (H) | Initial topography, corrected for total subsidence. |

The `PackVolume`, `PackSlice` and `PackColumn` type aliases can be used for dispatch, similar to `DataVolume`, `DataSlice` and `DataColumn`.

All visualization functions accept a `Pack`, as well as a `Header`/`Data` pair for backward compatibility. For example:

```julia
using CarboKitten.Output: load, load_slice
using CarboKitten.Visualization: sediment_profile

load("data/output/alcap-example.h5") do bundle
    sediment_profile(load_slice(bundle, :profile))
end
```

Visualization
-------------

The visualization of CarboKitten output is implemented in a [Julia package extension](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)). This is done so that `CarboKitten.jl` itself doesn't have to depend on `Makie.jl` (our main visualization tool), which has a large transient dependency stack. To make the `Visualization` extension of CarboKitten available, make sure to activate a Julia project where `Makie` is installed.

### Summary plot

The summary plot is a way to get an impression of your model output quickly. The `summary_plot` expects an HDF5 filename or a `Bundle` as input.

### Gallery

!!! note "TODO"
    add a gallery of available plotting routines here.

### Makie primer

`Makie.jl` is a visualization package that creates exceptionally good looking (publication quality) plots in both 2D and 3D. There are three back-ends for Makie:

- `CairoMakie` for publication quality vector graphics, writing to `SVG`, `PDF` or `PNG`.
- `GLMakie` has better run-time performance than `CairoMakie`, especially when dealing with larger datasets and/or 3D visualizations. However, `GLMakie` can only produce rasterized images, so `PNG`, `JPEG` or directly to screen for interactive use.
- `WGLMakie` for online publication using WebGL. If you want interactive plots, like 3D plots that you can rotate in the browser, this is the one to use. Fair warning: this is also the least stable back-end for Makie.

To work with Makie, you need to import one of the three back-end packages. In general, every plot available in Makie has two variants. One is a direct function for plotting:

```julia
using CairoMakie

x = randn(10)
y = randn(10)

scatter(x, y)
```

The other requires a bit more prep, but gives you more control.

```julia
fig = Figure()
ax = Axis(fig[1,1])
scatter!(ax, x, y)
```

Here, we create a figure explicitly, then create a new set of axes somewhere on the grid in the figure, and then plot on that set of axes. The plotting functions accepting an `Axis` argument actually modify an existing context, which is why these functions always end with an exclamation mark, in this case `scatter!`.

If you like to know more about Makie, their ["Getting started"](https://docs.makie.org/stable/tutorials/getting-started) is a good place to start.

### Coding standards

Implementations of visualizations should never use member access to obtain data, but use only the methods described above.
