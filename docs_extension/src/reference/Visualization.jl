"""
    facies_colormap(n_facies; facies_colors=nothing, include_nodeposit=true)

Construct a categorical colormap for discrete facies visualization.

# Arguments
- `n_facies`: Number of facies categories.
- `facies_colors`: Optional vector of colors. If `nothing`, a default palette is used.
- `include_nodeposit`: If `true`, prepends a white color representing facies ID = 0.

# Returns
- A Makie categorical colormap with `n_facies` entries (plus optional NoDeposit).

# Algorithm
1. Select base color palette (user-provided or default).
2. If `include_nodeposit == true`, prepend a white color.
3. Construct categorical colormap with consistent indexing.

# Units
Not applicable (categorical data).

# Notes
Color ordering must match facies indexing in input grids.
"""
function facies_colormap end

"""
    fence_plot(grid::AbstractArray{<:Integer,3}; kwargs...)

Generate a 3D fence diagram from a categorical voxel grid.

# Arguments
- `grid`: 3D array of integer facies IDs.
- `topk`: Optional 2D array giving the top-most valid index per column.
- `category_names`: Labels corresponding to facies IDs.
- `dx`: Horizontal spacing.
- `dz`: Vertical spacing.
- `colormap`: Optional Makie colormap.
- `facies_colors`: Optional explicit palette.
- `colorbar_label`: Label for the colorbar.
- `xs`, `ys`: Slice indices used to extract fence planes.
- `VE`: Vertical exaggeration factor.
- `fence_alpha`: Surface transparency.
- `shading`: Enable/disable lighting.
- `frame_lw`: Frame line width.
- `frame_color`: Frame color.

# Returns
- Makie `Figure` containing the rendered fence diagram.

# Algorithm
1. Apply `topk` mask to remove voxels above the preserved stratigraphic surface.
2. Replace zero-valued cells with `NaN` to prevent rendering.
3. Extract vertical slices at specified `xs` and `ys`.
4. Convert grid indices to physical coordinates using `dx`, `dz`, and `VE`.
5. Render each slice as a surface using categorical coloring.
6. Add frame outlines and categorical colorbar.

# Units
- `dx`, `dz` in meters.
- `VE` rescales vertical axis only for visualization.

# Notes
Grid indexing order must match spatial interpretation of axes.
"""
function fence_plot end

"""
    glamour_view!(ax::Makie.Axis3, header::Header, data::DataVolume; colormap=Reverse(:speed))

Render stacked stratigraphic surfaces into an existing 3D axis.

# Arguments
- `ax`: Target Makie 3D axis.
- `header`: Grid and model metadata.
- `data`: Data volume containing thickness or elevation fields.
- `colormap`: Colormap used for surface rendering.

# Returns
- `nothing`

# Algorithm
1. Extract thickness or elevation arrays from `data`.
2. Select representative surfaces (e.g. base, intermediate, top).
3. Convert grid coordinates to physical space.
4. Render surfaces using Makie `surface!`.
5. Add domain boundary lines.

# Units
Coordinates are typically in meters; horizontal axes may be scaled to kilometers for display.

# Notes
This method modifies the provided axis in place.
"""
function glamour_view! end

"""
    glamour_view!(ax::Makie.Axis3, header::Header, input, state; colormap=Reverse(:speed))

Render stacked stratigraphic surfaces from model state.

# Arguments
- `ax`: Target Makie axis.
- `header`: Model metadata.
- `input`: Model configuration.
- `state`: Dynamic simulation state.
- `colormap`: Colormap used for rendering.

# Returns
- `nothing`

# Algorithm
1. Reconstruct thickness fields from `state`.
2. Select representative timesteps.
3. Convert to physical coordinates.
4. Render surfaces using Makie.

# Units
Consistent with model configuration (typically meters).

# Notes
Used for visualization directly from simulation output.
"""
function glamour_view! end

"""
     glamour_view(header::Header, data::DataVolume; colormap=Reverse(:speed))

Create a new figure and render a glamour view from precomputed data.

# Arguments
- `header`: Model metadata.
- `data`: Data volume.
- `colormap`: Colormap used for rendering.

# Returns
- Makie `Figure`

# Algorithm
1. Create a new figure and 3D axis.
2. Call `glamour_view!` to render surfaces.

# Units
Inherited from input data.

# Notes
Wrapper function for figure creation.
"""
function glamour_view end

"""
    glamour_view(header::Header, input, state; colormap=Reverse(:speed))

Create a new figure and render a glamour view from model state.

# Arguments
- `header`: Model metadata.
- `input`: Model configuration.
- `state`: Simulation state.
- `colormap`: Colormap used for rendering.

# Returns
- Makie `Figure`

# Algorithm
1. Create a new figure and 3D axis.
2. Call `glamour_view!` to render reconstructed surfaces.

# Units
Inherited from model configuration.

# Notes
Convenience wrapper for direct visualization.
"""
function glamour_view end

"""
    map_view(grid::AbstractArray{<:Integer,3}; kwargs...)

Render a horizontal categorical slice from a 3D voxel grid.

# Arguments
- `grid`: 3D array of integer category or facies IDs.
- `z0`: Vertical index of the slice to extract.
- `category_names`: Labels corresponding to category IDs.
- `category_colors`: Optional vector of colors defining the categorical palette.
- `colorbar_label`: Label displayed on the colorbar.

# Returns
- Makie `Figure` containing the rendered slice.

# Algorithm
1. Extract the horizontal slice `grid[z0, :, :]`.
2. Replace zero or invalid category IDs with `NaN` to prevent rendering.
3. Construct a categorical colormap using `category_colors` or defaults.
4. Render the slice using `heatmap!`.
5. Add a categorical colorbar with matching labels.

# Units
Grid values are categorical identifiers and carry no physical units.

# Notes
Indexing order must match the spatial interpretation of the grid axes.
"""
function map_view end

"""
    production_curve!(ax, input::I; facies_colors=nothing, max_depth=500.0u"m", facies_names=nothing) where I <: AbstractInput

Plot facies production rate as a function of depth into an existing axis.

# Arguments
- `ax`: Target Makie axis.
- `input`: Model input containing facies definitions and production parameters.
- `facies_colors`: Optional color palette.
- `max_depth`: Maximum depth used to evaluate production curves.
- `facies_names`: Optional labels for legend entries.

# Returns
- `nothing`

# Algorithm
1. Construct a depth vector from surface to `max_depth`.
2. For each facies:
   - Evaluate its production function across depth.
3. Plot production rate versus depth for each facies.
4. Add legend and axis labels.

# Units
- Depth in meters.
- Production rate typically in m/Myr.

# Notes
The function modifies the provided axis in place.
"""
function production_curve! end

"""
    production_curve(input::I; facies_colors=nothing, max_depth=500.0u"m", facies_names=nothing) where I <: AbstractInput

Create a new figure and plot facies production curves.

# Arguments
- `input`: Model input containing facies definitions.
- `facies_colors`: Optional color palette.
- `max_depth`: Maximum depth used to evaluate production curves.
- `facies_names`: Optional legend labels.

# Returns
- Makie `Figure`

# Algorithm
1. Create a new figure and axis.
2. Call `production_curve!` to render curves.

# Units
Inherited from input definitions.

# Notes
Wrapper function for figure creation.
"""
function production_curve end

"""
    production_curve!(ax, g::HDF5.Group; max_depth=-50.0u"m", facies_colors=nothing, facies_names=nothing)

Plot facies production curves using data stored in an HDF5 group.

# Arguments
- `ax`: Target Makie axis.
- `g`: HDF5 group containing facies metadata and production parameters.
- `max_depth`: Maximum depth used to evaluate production curves.
- `facies_colors`: Optional color palette.
- `facies_names`: Optional legend labels.

# Returns
- `nothing`

# Algorithm
1. Extract facies parameters from HDF5 group.
2. Construct a depth vector.
3. Evaluate production functions for each facies.
4. Plot curves and add legend.

# Units
- Depth in meters.
- Production rate in m/Myr.

# Notes
Provides a reconstruction-based alternative to input-driven plotting.
"""
function production_curve! end

"""
    sediment_profile!(ax::Axis, header::Header, data::DataSlice;
                      show_unconformities=true,
                      show_coeval_lines=true,
                      show_sealevel=true,
                      facies_colors=nothing)

Render a sediment-profile plot showing dominant facies, optional coeval lines, unconformities, and sea level.

# Arguments
- `ax`: Target Makie axis.
- `header`: Model metadata containing grid geometry and sea-level history.
- `data`: Data slice containing facies production and deposition fields.
- `show_unconformities`: Boolean or control flag for plotting unconformity lines.
- `show_coeval_lines`: Controls coeval-line plotting. Accepted types:
    - `Bool`: use default spacing
    - `Tuple{Int,Int}`: major/minor spacing
    - `Vector{Time}`: specific times
    - `Vector{Int}`: specific indices
- `show_sealevel`: If `true`, plots final sea-level elevation.
- `facies_colors`: Optional palette for facies coloring.

# Returns
- `(ax, hm)`:
    - `ax`: Updated Makie axis
    - `hm`: Handle to the facies heatmap

# Algorithm
1. Compute dominant facies per column:
   - Use `argmax` over production/deposition dimension.
   - Mask columns with zero total deposition (`tot <= 0`) using `NaN`.
2. Construct categorical colormap via `facies_colormap`.
3. Plot dominant facies using `profile_plot!`.
4. Dispatch coeval-line plotting based on the type of `show_coeval_lines`.
5. Overlay unconformities if enabled using `plot_unconformities`.
6. Plot final sea level as a horizontal line if enabled.
7. Set axis title and return axis and plot handle.

# Units
- Vertical coordinates in meters.
- Horizontal coordinates typically converted to kilometers.
- Facies values are categorical identifiers.

# Notes
This function modifies the provided axis in place and combines multiple visualization layers.
"""
function sediment_profile! end

"""
    sediment_profile(header::Header, data::DataSlice; kwargs...)

Create a new figure and render a sediment-profile plot.

# Arguments
- `header`: Model metadata.
- `data`: Data slice.
- `kwargs...`: Passed to `sediment_profile!`.

# Returns
- `(fig, ax, hm)`:
    - `fig`: Makie figure
    - `ax`: Axis
    - `hm`: Heatmap handle

# Algorithm
1. Create a new figure and axis.
2. Call `sediment_profile!` to render all layers.
3. Return figure, axis, and plot handle.

# Units
Inherited from input data.

# Notes
Convenience wrapper for figure creation.
"""
function sediment_profile end

"""
    stratigraphic_column_layers!(
        ax::Makie.Axis,
        state;
        x_index::Int=1,
        y_index::Int=1,
        facies_colors=nothing
    )

Render the preserved compacted stratigraphic layer stack for a single grid column.

# Arguments
- `ax`: Target Makie axis.
- `state`: Model state containing the preserved compaction archive.
- `x_index`: Horizontal x-index of the selected column.
- `y_index`: Horizontal y-index of the selected column.
- `facies_colors`: Optional categorical color palette.

# Returns
- `nothing`

# Algorithm
1. Extract the compacted layer stack at `(x_index, y_index)` from the state archive.
2. Iterate through layers in stratigraphic order.
3. Map facies identifiers to colors using the provided or default palette.
4. Draw each layer as a horizontal segment spanning its thickness.

# Units
- Vertical coordinates in meters.
- Facies identifiers are categorical.

# Notes
Operates directly on preserved, compacted layers without reconstruction.
"""
function stratigraphic_column_layers! end

"""
    stratigraphic_column!(ax::Makie.Axis, header::Header, state;
                          x_index=1,
                          y_index=1,
                          write_interval=1,
                          facies_colors=nothing)

Render a stratigraphic column from in-memory model state.

# Arguments
- `ax`: Target Makie axis.
- `header`: Model metadata containing grid and time information.
- `state`: Model state containing deposition history.
- `x_index`: Horizontal x-index of the selected column.
- `y_index`: Horizontal y-index of the selected column.
- `write_interval`: Interval stride used when converting state history.
- `facies_colors`: Optional categorical color palette.

# Returns
- `nothing`

# Algorithm
1. Extract a column at `(x_index, y_index)` from the model state.
2. Convert the state history to interval-based `DataColumn` representation.
3. Delegate plotting to the `DataColumn`-based stratigraphic-column method.

# Units
- Time in Myr where applicable.
- Facies identifiers are categorical.

# Notes
Acts as a bridge between state-based and data-based visualization workflows.
"""
function stratigraphic_column! end



"""
    stratigraphic_column(header::Header, column::DataColumn, facies::Int)

Render stratigraphic column from exported data.

# Arguments
- `header`: Metadata.
- `column`: Data column.
- `facies`: Facies index.

# Returns
- Makie `Figure`

# Algorithm
1. Extract facies intervals.
2. Plot vertical column.

# Units
Distance in meters.
"""
function stratigraphic_column end


"""
    summary_plot(fid::HDF5.File;
                 wheeler_smooth=(1,1),
                 show_unconformities=true,
                 facies_colors=nothing,
                 facies_names=nothing)

Create a multi-panel summary visualization from an open HDF5 file.

# Arguments
- `fid`: Open HDF5 file handle.
- `wheeler_smooth`: Smoothing kernel size for Wheeler plots.
- `show_unconformities`: Toggle for unconformity overlays.
- `facies_colors`: Optional categorical palette.
- `facies_names`: Optional legend labels.

# Returns
- Makie `Figure`

# Algorithm
1. Load header metadata and available datasets.
2. Construct component plots:
   - sediment profile
   - Wheeler diagram
   - sea level
   - topography
   - production curves
3. Arrange subplots into a single figure layout.

# Units
Inherited from stored simulation outputs.

# Notes
Automatically adapts to available datasets.
"""
function summary_plot end

"""
    summary_plot_from_state(header, state, input;
                            wheeler_smooth=(1,1),
                            show_unconformities=true,
                            facies_colors=nothing,
                            facies_names=nothing)

Create a multi-panel summary visualization directly from in-memory model state.

# Arguments
- `header`: Model metadata.
- `state`: Model state.
- `input`: Model configuration.
- `wheeler_smooth`: Smoothing kernel size.
- `show_unconformities`: Toggle for unconformity overlays.
- `facies_colors`: Optional categorical palette.
- `facies_names`: Optional legend labels.

# Returns
- Makie `Figure`

# Algorithm
1. Convert state to interval-based `DataSlice` using exact reconstruction.
2. Build component plots directly from state and input.
3. Assemble plots into a multi-panel layout.

# Units
- Time in Myr.
- Distance in meters.
- Facies identifiers are categorical.

# Notes
Avoids intermediate file export.
"""
function summary_plot_from_state end

"""
    wheeler_column(state; kwargs...)

Construct a categorical Wheeler-style time column from state history.

# Arguments
- `state`: Model state.
- `x`: Horizontal x-index.
- `y`: Optional y-index.
- `Δi`, `Δj`: Spatial averaging window size.
- `category_names`: Labels for categories.
- `category_colors`: Categorical palette.
- `selector`: Function mapping local conditions to category ID.
- `nondep_eps`: Threshold for non-deposition.
- `title_str`: Plot title.
- `colorbar_label`: Colorbar label.

# Returns
- Makie `Figure`

# Algorithm
1. Extract spatial window defined by `(x, y, Δi, Δj)`.
2. Aggregate facies strength, water depth, and energy per timestep.
3. Apply `selector` to assign category IDs.
4. Render time-column as categorical heatmap.

# Units
Categorical values only.

# Notes
Selector controls classification logic.
"""
function wheeler_column end

"""
    sediment_accumulation!(ax::Axis, header::Header, data::DataSlice;
                           smooth_size=(3,11),
                           colormap=Reverse(:curl),
                           range=(-100.0u"m/Myr", 100.0u"m/Myr"))

Plot interval-based sediment accumulation as a Wheeler heatmap.

# Arguments
- `ax`: Target axis.
- `header`: Model metadata.
- `data`: Data slice.
- `smooth_size`: Smoothing kernel size.
- `colormap`: Colormap.
- `range`: Value range for visualization.

# Returns
- `nothing`

# Algorithm
1. Compute net accumulation (deposition − erosion).
2. Normalize by interval duration.
3. Apply smoothing filter.
4. Render as heatmap over space and time.

# Units
- Accumulation in m/Myr.
- Distance in meters or kilometers.

# Notes
Highlights depositional vs erosional phases.
"""
function sediment_accumulation! end

"""
    dominant_facies!(ax::Axis, header::Header, data::DataSlice;
                     smooth_size=(3,11),
                     facies_colors=nothing,
                     nondep_eps=0.0u"m")

Plot dominant facies through time as a Wheeler heatmap.

# Arguments
- `ax`: Target axis.
- `header`: Model metadata.
- `data`: Data slice.
- `smooth_size`: Smoothing kernel size.
- `facies_colors`: Optional palette.
- `nondep_eps`: Threshold for non-deposition.

# Returns
- `nothing`

# Algorithm
1. Determine dominant facies per interval.
2. Mask intervals below `nondep_eps`.
3. Apply smoothing.
4. Render categorical heatmap.

# Units
Categorical values.

# Notes
Represents temporal facies evolution.
"""
function dominant_facies! end

"""
    wheeler_diagram!(ax1::Axis, ax2::Axis, header::Header, data::DataSlice; kwargs...)

Render Wheeler-style summary plots into existing axes.

# Arguments
- `ax1`, `ax2`: Target axes.
- `header`: Model metadata.
- `data`: Data slice.
- `kwargs...`: Passed to subroutines.

# Returns
- `nothing`

# Algorithm
1. Call `sediment_accumulation!` for accumulation field.
2. Call `dominant_facies!` for categorical field.
3. Add shared colorbars.

# Units
- Time in Myr.
- Distance in meters or kilometers.

# Notes
Low-level rendering routine.
"""
function wheeler_diagram! end

"""
    wheeler_diagram(header::Header, data::DataSlice; kwargs...)

Create a new figure containing Wheeler-style visualizations.

# Arguments
- `header`: Model metadata.
- `data`: Data slice.
- `kwargs...`: Passed to `wheeler_diagram!`.

# Returns
- Makie `Figure`

# Algorithm
1. Create figure and axes.
2. Call `wheeler_diagram!`.
3. Add colorbars and layout adjustments.

# Units
Inherited from data.

# Notes
High-level wrapper.
"""
function wheeler_diagram end
