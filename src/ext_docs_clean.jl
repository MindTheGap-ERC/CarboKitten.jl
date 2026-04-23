module ExtDocs

using CarboKitten

# --- stubs (generic functions / types) ---

function build_vertical_coordinates end
function facies_colormap end
function fence_plot end

function glamour_view! end
function glamour_view end

function build_selected_surfaces end
function map_view end

function production_curve! end
function production_curve end

function resample_to_regular_grid end
function elevation end
function explode_quad_vertices end
function plot_unconformities end

function coeval_lines! end

function sediment_profile! end
function sediment_profile end

function elevation_physical end

# qualified name present in @docs (define in that module)
function stratigraphic_column_layers! end

function stratigraphic_column! end

function summary_plot end
function summary_plot_from_state end
function wheeler_column end

function sediment_accumulation! end
function dominant_facies! end

function wheeler_diagram! end
function wheeler_diagram end

function init end
function run_model end

function box_axes end
function n_steps end
function set_attribute end

function cementation_factor end
function adaptive_transporter end
function disintegrator end
function transporter end

function write_header end
function production end
function rules end
function step_ca end
function get_logger end

struct CompactionState end

function insolation end

function depth_multiplier end
function time_multiplier end
function porosity_at_depth end
function compacted_thickness end

function compact_column! end
function update_sediment_height_from_compaction! end
function merge_layers! end

function extract_compacted_layer_maps! end
function project_compacted_layers_to_cube! end
function bury_deposition! end

function capped_production end
function uniform_production end

function initial_state end
function initial_topography end
function water_depth end

function facies_props_depthbins end
function header_from_input_exact end
function datacolumn_from_state end
function dataslice_from_state_exact end
function datavolume_from_state end

function data_kind end
function group_datasets end
function data_header end
function read_header end
function read_data end

function accumulate end  # for docs listing

function stratigraphic_column end

function data_export end

function extract_sac end
function extract_sc end
function extract_wd end

function step! end

function parse_slice end

function new_output end
function add_data_set end

function write_sediment_thickness end
function write_production end
function write_disintegration end
function write_deposition end

function state_writer end
function frame_writer end

function facies_percent_from_ids end
function classify_block end
function env_selector end
function build_environment_grid end
function encode_environments end

function write_environment_block! end
function write_layer_archive! end

function make_header end

struct H5Output end

function patch_stats_3D_from_ids end
function run_once end
function run_single_to_csv end

function push_sediment! end
function pop_sediment! end
function peek_sediment end

function pde_stencil end
function advection_coef! end
function max_dt end
function transport_dC! end
function transport! end
function transport end

function _wave_velocity_at_depth end
function wavelength_from_period end
function wave_physics_at_cell end

# also needed for qualified entry
function run_model end

end





@doc """
build_vertical_coordinates(thickness_layers)

Reconstruct cumulative layer-top and layer-base elevation surfaces from a stack of layer thickness rasters.

# Arguments
- `thickness_layers`: Ordered collection of 2D arrays representing layer thickness. The first element corresponds to the oldest (deepest) layer and the last to the youngest.

# Returns
- `(z_top, z_base)`:
    - `z_top[i]`: Elevation of the top of layer `i`
    - `z_base[i]`: Elevation of the base of layer `i`

# Algorithm
1. Initialize a cumulative elevation field `z = 0` for all horizontal grid cells.
2. Iterate through layers in stratigraphic order:
   - Assign current `z` as the base of the layer.
   - Compute top as `z + thickness`.
   - Update `z` to the new top elevation.
3. Store base and top surfaces for each layer.

# Units
All thicknesses and elevations are assumed to be in meters.

# Notes
The function assumes spatial alignment across all input layers.
""" -> build_vertical_coordinates

@doc """
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
""" -> facies_colormap

@doc """
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
""" -> fence_plot

@doc """
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
""" -> glamour_view!

@doc """
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
""" -> glamour_view!

@doc """
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
""" -> glamour_view

@doc """
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
""" -> glamour_view

@doc """
build_selected_surfaces(state, input, selected_steps)

Reconstruct compaction-adjusted palaeo-surfaces for selected timesteps.

# Arguments
- `state`: Model state containing preserved layer archive.
- `input`: Model configuration including compaction relationships.
- `selected_steps`: Vector of timestep indices.

# Returns
- Array of reconstructed surface elevations, one per timestep.

# Algorithm
1. For each timestep and grid column:
   - Select layers deposited up to that timestep.
2. Compute compacted thickness using burial depth and facies-specific compaction laws.
3. Sum compacted thickness vertically to obtain surface elevation.

# Units
All distances in meters.

# Notes
Accounts for burial compaction effects.
""" -> build_selected_surfaces

@doc """
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
""" -> map_view

@doc """
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
""" -> production_curve!

@doc """
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
""" -> production_curve

@doc """
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
""" -> production_curve!

@doc """
resample_to_regular_grid(facies_layers, thickness_layers, dz)

Convert irregular preserved layers into a regular vertical voxel grid.

# Arguments
- `facies_layers`: Vector of 2D arrays containing facies IDs.
- `thickness_layers`: Vector of 2D arrays containing layer thickness.
- `dz`: Target vertical resolution.

# Returns
- 3D categorical voxel grid.

# Algorithm
1. Compute total column thickness to determine required voxel count.
2. For each grid column:
   - Track cumulative depth.
   - Map each layer interval to voxel indices using `dz`.
3. Assign facies IDs to corresponding voxels.

# Units
- Thickness and `dz` in meters.
- Output grid is unitless categorical.

# Notes
Assumes layers are vertically ordered and spatially aligned.
""" -> resample_to_regular_grid

@doc """
thickness_field(d::DataSlice)

Select the appropriate thickness field for profile plotting.

# Arguments
- `d`: Data slice containing thickness-related fields.

# Returns
- 2D or 3D thickness array.

# Algorithm
1. Check for compacted thickness field.
2. If present, return it.
3. Otherwise return `sediment_thickness`.

# Units
Thickness in meters.

# Notes
Prefers physically corrected (compacted) thickness when available.
\"\"\"
thickness_field(d::DataSlice) =



\"\"\"
    aligned_time(header::Header, data::DataSlice)

Return a time axis consistent with interval-based data.

# Arguments
- `header`: Model metadata containing time axis.
- `data`: Data slice containing interval-based thickness.

# Returns
- Time vector aligned with data intervals.

# Algorithm
1. Determine number of intervals from thickness field.
2. Return first `Nt-1` values of the header time vector.

# Units
Time typically in Myr.

# Notes
Ensures consistency between exported interval data and time coordinates.
\"\"\"
aligned_time(header::Header, data::DataSlice) = begin



\"\"\"
    elevation(h::Header, d::DataSlice)

Reconstruct elevation history for sediment profile visualization.

# Arguments
- `h`: Model metadata including initial conditions.
- `d`: Data slice containing thickness history.

# Returns
- Elevation array over time.

# Algorithm
1. Extract thickness field.
2. Combine:
   - initial topography
   - accumulated sediment thickness
   - subsidence adjustments
3. Compute elevation through time.

# Units
Elevation in meters.

# Notes
Used for section-based stratigraphic plots.
""" -> elevation

@doc """
explode_quad_vertices(v::Array{Float64,3})

Convert structured quad grid into explicit mesh representation.

# Arguments
- `v`: Structured grid of vertex coordinates.

# Returns
- Tuple of `(vertices, faces)` suitable for mesh plotting.

# Algorithm
1. Expand grid into explicit vertex list.
2. Construct face connectivity for quadrilateral elements.

# Units
Inherited from input coordinates.

# Notes
Required for Makie mesh rendering.
""" -> explode_quad_vertices

@doc """
plot_unconformities(ax::Axis, header::Header, data::DataSlice, ::Nothing; kwargs...)

Overlay unconformity lines on a sediment profile plot.

# Arguments
- `ax`: Target Makie axis.
- `header`: Model metadata.
- `data`: Data slice.
- `::Nothing`: Placeholder argument.
- `kwargs...`: Plot styling options.

# Returns
- `nothing`

# Algorithm
1. Compute water depth from sea level and elevation.
2. Identify subaerial intervals.
3. Skeletonize binary mask into line segments.
4. Plot lines on axis.

# Units
- Distance in meters.
- Horizontal axis may be converted to kilometers.

# Notes
Used to highlight sequence boundaries and exposure surfaces.
""" -> plot_unconformities

@doc """
coeval_lines!(ax, header, data, flag::Bool)

Plot coeval lines using a default timestep specification.

# Arguments
- `ax`: Target Makie axis.
- `header`: Model metadata.
- `data`: Data slice.
- `flag`: If true, triggers plotting.

# Returns
- `nothing`

# Algorithm
1. If `flag` is true, delegate to tuple-based overload.

# Units
Inherited from data.

# Notes
Convenience wrapper.
""" -> coeval_lines!

@doc """
coeval_lines!(ax, header, data, ts::Vector{Time})

Plot coeval lines for specified times.

# Arguments
- `ax`: Target axis.
- `header`: Metadata.
- `data`: Data slice.
- `ts`: Vector of times.

# Returns
- `nothing`

# Algorithm
1. Convert times to indices using `searchsortedfirst`.
2. Clamp indices.
3. Delegate to index-based overload.

# Units
Time in Myr.

# Notes
Handles non-uniform time spacing.
""" -> coeval_lines!

@doc """
coeval_lines!(ax, header, data, ticset::Tuple{Int,Int})

Plot major and minor coeval lines.

# Arguments
- `ticset`: Tuple specifying number of major and minor intervals.

# Returns
- `nothing`

# Algorithm
1. Generate evenly spaced indices.
2. Separate into major and minor sets.
3. Plot with different styles.

# Notes
Used for visual grid structure.
""" -> coeval_lines!

@doc """
coeval_lines!(ax, header, data, idx::Vector{Int})

Plot coeval lines at specified indices.

# Arguments
- `idx`: Indices of timesteps.

# Returns
- `nothing`

# Algorithm
1. Extract x-coordinates and elevation.
2. Clamp indices.
3. Plot elevation profiles for each index.

# Units
- Horizontal in km.
- Vertical in meters.

# Notes
Lowest-level implementation.
""" -> coeval_lines!

@doc """
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
""" -> sediment_profile!

@doc """
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
""" -> sediment_profile

@doc """
elevation_physical(h::Header, d::DataSlice)

Compute physical elevation of the sediment profile through time.

# Arguments
- `h`: Model metadata containing initial topography.
- `d`: Data slice containing thickness and subsidence history.

# Returns
- `el`: Array of elevation values through time.

# Algorithm
1. Extract thickness field using `thickness_field`.
2. Retrieve initial topography for the selected slice.
3. Extract cumulative subsidence history.
4. Compute elevation as:
   `elevation = initial_topography + thickness - subsidence`
5. Return elevation array.

# Units
- Elevation and thickness in meters.
- Subsidence in meters.

# Notes
Represents elevation in physical (not stratigraphic) space, accounting for subsidence.
""" -> elevation_physical

@doc """
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
""" -> stratigraphic_column_layers!

@doc """
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
""" -> stratigraphic_column!

@doc """
summary_plot(filename::AbstractString; kwargs...)

Create a multi-panel summary figure from an HDF5 simulation output file.

# Arguments
- `filename`: Path to HDF5 file.
- `kwargs...`: Passed to the file-handle-based method.

# Returns
- Makie `Figure`

# Algorithm
1. Open the HDF5 file in read mode.
2. Delegate processing to `summary_plot(fid; kwargs...)`.

# Units
Inherited from stored datasets.

# Notes
Thin wrapper for file-based workflow.
\"\"\"
summary_plot(filename::AbstractString; kwargs...) = h5open(fid->summary_plot(fid; kwargs...), filename, "r")



\"\"\"
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
""" -> summary_plot

@doc """
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
""" -> summary_plot_from_state

@doc """
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
""" -> wheeler_column

@doc """
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
""" -> sediment_accumulation!

@doc """
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
""" -> dominant_facies!

@doc """
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
""" -> wheeler_diagram!

@doc """
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
""" -> wheeler_diagram

@doc """
init()

Initialize the CarboKitten runtime environment.

# Arguments
- None

# Returns
- `nothing`

# Algorithm
1. Configure global logging using `global_logger` and `TerminalLogger`.
2. Prepare runtime settings required before model execution.

# Units
Not applicable.

# Notes
Executed once at startup before any model run.
""" -> init

@doc """
run_model

Entry point for executing a CarboKitten simulation.

# Arguments
- None (dispatched on input types elsewhere)

# Returns
- Simulation result defined by dispatched implementation.

# Algorithm
1. Dispatch to method implementations based on input type.
2. Execute the full model workflow (initialization, stepping, output).

# Units
- Time typically in Myr.
- Distance in meters.

# Notes
This is a generic function; behavior is defined by method specialization.
""" -> run_model

@doc """
box_axes(box::Box)

Return the spatial axes associated with a simulation domain.

# Arguments
- `box`: Domain definition containing grid geometry.

# Returns
- Tuple or structure describing spatial axes.

# Algorithm
1. Extract axis definitions from the `Box` structure.
2. Return coordinate arrays for each dimension.

# Units
Coordinates typically in meters.

# Notes
Used to define spatial grids for model components.
""" -> box_axes

@doc """
n_writes(time::TimeProperties)

Return the number of output write intervals.

# Arguments
- `time`: Time configuration.

# Returns
- Integer number of write intervals.

# Algorithm
1. Return `time.steps`.

# Units
Dimensionless count.

# Notes
Defines how many output snapshots are produced.
\"\"\"
n_writes(time::TimeProperties) = time.steps



\"\"\"
    n_steps

Return the total number of simulation timesteps.

# Arguments
- None (dispatched elsewhere)

# Returns
- Integer number of timesteps.

# Algorithm
1. Dispatch to implementation using time configuration.
2. Typically derived from `time_axis` or `n_writes`.

# Units
Dimensionless count.

# Notes
Generic function; implementation defined elsewhere.
""" -> n_steps

@doc """
time_axis(time::TimeProperties)

Construct the time axis corresponding to model output intervals.

# Arguments
- `time`: Time configuration.

# Returns
- Range of time values.

# Algorithm
1. Compute `(0:n_writes(time)) * time.Δt + time.t0`.

# Units
Time in Myr.

# Notes
Includes both initial and final times.
\"\"\"
time_axis(time::TimeProperties) = (0:n_writes(time)) .* time.Δt .+ time.t0



\"\"\"
    set_attribute

Assign metadata attributes to output objects.

# Arguments
- None (dispatched elsewhere)

# Returns
- `nothing`

# Algorithm
1. Dispatch to output-backend-specific implementations.

# Units
Depends on attribute type.

# Notes
Generic function used during output writing.
""" -> set_attribute

@doc """
courant_max(::Type{Val{:RK4}})

Return the Courant stability limit for RK4 transport integration.

# Arguments
- `Val{:RK4}`: Integration scheme identifier.

# Returns
- `2.0`

# Algorithm
1. Return fixed stability coefficient for RK4.

# Units
Dimensionless.

# Notes
Used to constrain timestep selection.
\"\"\"
courant_max(::Type{Val{:RK4}}) = 2.0



\"\"\"
    courant_max(::Type{Val{:forward_euler}})

Return the Courant stability limit for forward Euler integration.

# Arguments
- `Val{:forward_euler}`: Integration scheme identifier.

# Returns
- `1.0`

# Algorithm
1. Return fixed stability coefficient for forward Euler.

# Units
Dimensionless.

# Notes
More restrictive than RK4.
\"\"\"
courant_max(::Type{Val{:forward_euler}}) = 1.0



\"\"\"
    transport_solver(f, _)

Return the provided transport function unchanged.

# Arguments
- `f`: Transport function.
- `_`: Ignored argument.

# Returns
- `f`

# Algorithm
1. Identity mapping.

# Units
Not applicable.

# Notes
Fallback solver configuration.
\"\"\"
transport_solver(f, _) = f



\"\"\"
    transport_solver(::Type{Val{:RK4}}, box)

Construct an RK4 transport solver for the given domain.

# Arguments
- `Val{:RK4}`: Solver type.
- `box`: Domain definition.

# Returns
- RK4 solver function.

# Algorithm
1. Call `runge_kutta_4` with unit type and domain.

# Units
Distance units derived from `box`.

# Notes
Uses unit-aware solver construction.
\"\"\"
transport_solver(::Type{Val{:RK4}}, box) = runge_kutta_4(typeof(1.0u"m"), box)



\"\"\"
    transport_solver(::Type{Val{:forward_euler}}, _)

Return the forward Euler transport solver.

# Arguments
- `Val{:forward_euler}`: Solver type.
- `_`: Ignored.

# Returns
- `forward_euler`

# Algorithm
1. Return solver function reference.

# Units
Not applicable.

# Notes
Simplest integration scheme.
\"\"\"
transport_solver(::Type{Val{:forward_euler}}, _) = forward_euler



\"\"\"
    cementation_factor(input::AbstractInput)

Compute facies-dependent cementation scaling factor.

# Arguments
- `input`: Model input containing facies parameters.

# Returns
- Array or scalar cementation factor.

# Algorithm
1. Evaluate exponential/logarithmic relationships based on input parameters.
2. Return scaling factor used in transport or deposition calculations.

# Units
Dimensionless.

# Notes
Controls lithification effects.
""" -> cementation_factor

@doc """
adaptive_transporter(input)

Construct an adaptive sediment transport operator.

# Arguments
- `input`: Model configuration.

# Returns
- Transport operator function.

# Algorithm
1. Select transport solver via `transport_solver`.
2. Compute stability limits using `courant_max`.
3. Evaluate transport rates using water depth and forcing.
4. Adapt timestep or fluxes based on stability criteria.

# Units
- Distance in meters.
- Time in Myr.

# Notes
Core transport routine for active-layer dynamics.
""" -> adaptive_transporter

@doc """
disintegrator(input)

Compute sediment disintegration or erosion rates.

# Arguments
- `input`: Model configuration.

# Returns
- Disintegration rate field.

# Algorithm
1. Compute water depth.
2. Evaluate erosion rates using facies-dependent parameters.
3. Apply limiting conditions (e.g. minimum thickness).

# Units
Rates typically in m/Myr.

# Notes
Represents breakdown of consolidated material.
""" -> disintegrator

@doc """
transporter(input)

Construct the sediment transport operator.

# Arguments
- `input`: Model configuration.

# Returns
- Transport operator.

# Algorithm
1. Delegate to `adaptive_transporter`.

# Units
Inherited from transport model.

# Notes
High-level wrapper.
""" -> transporter

@doc """
write_header(input::AbstractInput, output::AbstractOutput)

Write simulation metadata to the output backend.

# Arguments
- `input`: Model configuration.
- `output`: Output backend.

# Returns
- `nothing`

# Algorithm
1. Extract metadata from input.
2. Convert units where required.
3. Write attributes using `set_attribute`.

# Units
- Distance in meters.
- Time in Myr.

# Notes
Executed once per simulation.
""" -> write_header

@doc """
production(input::AbstractInput)

Compute carbonate production using the CA-controlled production pathway.

# Arguments
- `input`: Model configuration containing facies definitions and production controls.

# Returns
- Production rate field.

# Algorithm
1. Evaluate production using depth-dependent and time-dependent modifiers.
2. Apply CA-based constraints where applicable.
3. Combine facies contributions into a production field.

# Units
Production rates typically in m/Myr.

# Notes
Uses updated depth/time production formulation instead of fixed attenuation laws.
""" -> production

@doc """
rules(BT, facies, ca_priority, ca, i::CartesianIndex)

Evaluate activation and viability rules for a cellular automaton update.

# Arguments
- `BT`: Boundary type.
- `facies`: Current facies field.
- `ca_priority`: Priority ordering of CA rules.
- `ca`: Cellular automaton configuration.
- `i`: Target cell index.

# Returns
- Boolean or categorical state indicating rule activation.

# Algorithm
1. Extract neighbourhood around index `i`.
2. Apply facies-specific neighbourhood radius.
3. Evaluate rule priority ordering.
4. Return activation result.

# Units
Categorical identifiers.

# Notes
Encapsulates reusable CA rule evaluation logic.
""" -> rules

@doc """
step_ca(box::Box{BT}, facies) where {BT<:Boundary{2}}

Advance the cellular automaton by one timestep.

# Arguments
- `box`: Spatial domain and boundary conditions.
- `facies`: Current facies field.

# Returns
- Updated facies field.

# Algorithm
1. Iterate over all grid cells.
2. For each cell:
   - Evaluate neighbourhood using `rules`.
3. Update facies using CA transition logic.
4. Apply boundary conditions via `box`.

# Units
Categorical values.

# Notes
Core CA update step for facies evolution.
""" -> step_ca

@doc """
get_logger(input::AbstractInput)

Construct a logger for model execution.

# Arguments
- `input`: Model configuration containing logging settings.

# Returns
- Logger instance.

# Algorithm
1. Open output streams if required.
2. Configure logging levels using `MinLevelLogger`.
3. Combine file and terminal loggers.

# Units
Not applicable.

# Notes
Controls runtime diagnostics and output logging.
""" -> get_logger

@doc """
n_facies(input::AbstractInput)

Return the number of facies defined in the model.

# Arguments
- `input`: Model configuration.

# Returns
- Integer number of facies.

# Algorithm
1. Return `length(input.facies)`.

# Units
Dimensionless.

# Notes
Used throughout the model to size arrays and loops.
\"\"\"
n_facies(input::AbstractInput) = length(input.facies)



\"\"\"
    write_header(input::AbstractInput, output::AbstractOutput)

Write facies-related metadata to the output backend.

# Arguments
- `input`: Model configuration.
- `output`: Output backend.

# Returns
- `nothing`

# Algorithm
1. Iterate over facies definitions.
2. Write attributes using `set_attribute`.
3. Store facies metadata in output.

# Units
Categorical values.

# Notes
Facies definitions are written once at initialization.
""" -> write_header

@doc """
CompactionState(nx::Int, ny::Int, nzmax::Int)

Allocate state for compaction tracking.

# Arguments
- `nx`, `ny`: Horizontal grid dimensions.
- `nzmax`: Maximum number of preserved layers.

# Returns
- Compaction state object.

# Algorithm
1. Allocate arrays for:
   - layer thickness
   - facies IDs
   - compaction variables
2. Initialize storage with appropriate integer types.

# Units
Distances in meters.

# Notes
Stores preserved stratigraphy for compaction calculations.
""" -> CompactionState

@doc """
insolation(input::AbstractInput)

Compute insolation forcing.

# Arguments
- `input`: Model configuration.

# Returns
- Insolation field.

# Algorithm
1. Evaluate insolation function from input parameters.
2. Return spatial or temporal insolation values.

# Units
Depends on forcing definition.

# Notes
Used as a driver for production.
""" -> insolation

@doc """
write_header(input::AbstractInput, output::AbstractOutput)

Write production-related metadata to output.

# Arguments
- `input`: Model configuration.
- `output`: Output backend.

# Returns
- `nothing`

# Algorithm
1. Write production parameters.
2. Store insolation forcing.
3. Record time axis.

# Units
- Distance in meters.
- Time in Myr.

# Notes
Executed during initialization.
""" -> write_header

@doc """
depth_multiplier(f, w)

Evaluate depth-dependent production modifier.

# Arguments
- `f`: Facies parameters.
- `w`: Water depth.

# Returns
- Scaling factor.

# Algorithm
1. Convert depth to numeric form if unitful.
2. Evaluate depth response curve.
3. Return multiplier.

# Units
Depth in meters.

# Notes
Controls depth limitation of carbonate production.
""" -> depth_multiplier

@doc """
time_multiplier(f, t)

Evaluate time-dependent production modifier.

# Arguments
- `f`: Facies parameters.
- `t`: Time.

# Returns
- Scaling factor.

# Algorithm
1. Evaluate temporal modulation function.
2. Return multiplier.

# Units
Time in Myr.

# Notes
Used for time-varying production scenarios.
""" -> time_multiplier

@doc """
porosity_at_depth(f, z)

Compute porosity as a function of burial depth.

# Arguments
- `f`: Facies parameters.
- `z`: Depth.

# Returns
- Porosity value.

# Algorithm
1. Evaluate compaction curve `f.compaction_curve(z)`.
2. Validate result.
3. Return porosity.

# Units
Depth in meters.

# Notes
Defines compaction behaviour.
""" -> porosity_at_depth

@doc """
compacted_thickness(H0, ϕ0, ϕ)

Compute compacted thickness from initial thickness and porosity change.

# Arguments
- `H0`: Initial thickness.
- `ϕ0`: Initial porosity.
- `ϕ`: Final porosity.

# Returns
- Compacted thickness.

# Algorithm
1. Compute thickness reduction from porosity change.
2. Enforce non-negative thickness.

# Units
Distance in meters.

# Notes
Core compaction relation.
""" -> compacted_thickness

@doc """
compact_column!(input, state)

Apply compaction to all layers in each column.

# Arguments
- `input`: Model configuration.
- `state`: Model state.

# Returns
- `nothing`

# Algorithm
1. Loop over columns.
2. Traverse layers from top to bottom.
3. Update thickness using `porosity_at_depth` and `compacted_thickness`.

# Units
Distance in meters.

# Notes
Updates state in place.
""" -> compact_column!

@doc """
update_sediment_height_from_compaction!(state)

Update surface elevation after compaction.

# Arguments
- `state`: Model state.

# Returns
- `nothing`

# Algorithm
1. Recompute cumulative thickness.
2. Update sediment height field.

# Units
Distance in meters.

# Notes
Maintains consistency between compaction and elevation.
""" -> update_sediment_height_from_compaction!

@doc """
merge_layers!(state)

Merge thin or redundant layers in the compaction archive.

# Arguments
- `state`: Model state.

# Returns
- `nothing`

# Algorithm
1. Identify merge candidates.
2. Combine thickness and facies.
3. Update layer indices.

# Units
Distance in meters.

# Notes
Reduces archive size.
""" -> merge_layers!

@doc """
extract_compacted_layer_maps!(state, wdepth_hist, energy_hist)

Extract layer-wise maps from compacted archive.

# Arguments
- `state`: Model state.
- `wdepth_hist`: Water-depth history.
- `energy_hist`: Energy history.

# Returns
- `nothing`

# Algorithm
1. Loop over layers.
2. Extract thickness and facies.
3. Associate environmental fields.
4. Store maps.

# Units
- Depth in meters.
- Energy in model units.

# Notes
Prepares data for export or visualization.
""" -> extract_compacted_layer_maps!

@doc """
project_compacted_layers_to_cube!(state, Δz, wdepth_hist, energy_hist, block_cube, block_wdepth, block_energy, block_topk)

Project compacted layers onto a regular voxel grid.

# Arguments
- `state`: Model state.
- `Δz`: Vertical resolution.
- `wdepth_hist`, `energy_hist`: Environmental histories.
- `block_cube`, `block_wdepth`, `block_energy`, `block_topk`: Output arrays.

# Returns
- `nothing`

# Algorithm
1. Loop over columns and layers.
2. Convert thickness to voxel indices.
3. Fill voxel arrays with facies and environmental values.

# Units
Distance in meters.

# Notes
Used for 3D export.
""" -> project_compacted_layers_to_cube!

@doc """
bury_deposition!(state, deposited, step)

Insert newly deposited material into compaction archive.

# Arguments
- `state`: Model state.
- `deposited`: Deposition field.
- `step`: Current timestep.

# Returns
- `nothing`

# Algorithm
1. Append new layer to archive.
2. Store thickness and facies.
3. Validate consistency.

# Units
Distance in meters.

# Notes
Executed each timestep.
""" -> bury_deposition!

@doc """
production_rate(f, w, t)

Compute facies-specific production rate.

# Arguments
- `f`: Facies parameters.
- `w`: Water depth.
- `t`: Time.

# Returns
- Production rate.

# Algorithm
1. Evaluate depth multiplier.
2. Evaluate time multiplier.
3. Combine factors.

# Units
m/Myr.

# Notes
Updated production formulation.
\"\"\"
production_rate(f, w, t) =



\"\"\"
    capped_production(f, w, t, dt)

Limit production rate over timestep.

# Arguments
- `f`: Facies parameters.
- `w`: Water depth.
- `t`: Time.
- `dt`: Timestep.

# Returns
- Limited production.

# Algorithm
1. Compute production rate.
2. Clamp using timestep constraints.

# Units
Distance in meters.

# Notes
Ensures numerical stability.
""" -> capped_production

@doc """
uniform_production(input::AbstractInput)

Compute spatially uniform production.

# Arguments
- `input`: Model configuration.

# Returns
- Production field.

# Algorithm
1. Evaluate production at all grid points.
2. Apply uniform forcing.

# Units
m/Myr.

# Notes
Baseline production mode.
""" -> uniform_production

@doc """
initial_state(input::AbstractInput)

Construct the initial model state.

# Arguments
- `input`: Model configuration containing grid, forcing, and initial conditions.

# Returns
- State object.

# Algorithm
1. Allocate state arrays (e.g. elevation, thickness, facies fields).
2. Initialize fields using input configuration.
3. Return constructed `State`.

# Units
- Distance in meters.

# Notes
Defines the starting condition for the simulation.
""" -> initial_state

@doc """
initial_topography(input::AbstractInput)

Compute the initial topographic surface.

# Arguments
- `input`: Model configuration.

# Returns
- Elevation field.

# Algorithm
1. Evaluate initial topography definition from input.
2. Return elevation field over the model grid.

# Units
Elevation in meters.

# Notes
Represents pre-depositional surface.
""" -> initial_topography

@doc """
water_depth(input::AbstractInput)

Compute water depth through time.

# Arguments
- `input`: Model configuration containing sea level and topography.

# Returns
- Water-depth field.

# Algorithm
1. Retrieve initial topography.
2. Evaluate sea level as a function of time.
3. Compute water depth as:
   `water_depth = sea_level - topography`
4. Return water-depth field.

# Units
Distance in meters.

# Notes
Positive values indicate submerged conditions.
""" -> water_depth

@doc """
write_header(input::AbstractInput, output::AbstractOutput)

Write water-depth-related metadata to the output backend.

# Arguments
- `input`: Model configuration.
- `output`: Output backend.

# Returns
- `nothing`

# Algorithm
1. Extract spatial axes using `box_axes`.
2. Compute time axis using `time_axis`.
3. Retrieve initial topography.
4. Write attributes using `set_attribute`.

# Units
- Distance in meters.
- Time in Myr.

# Notes
Stores spatial and temporal reference information for output data.
""" -> write_header

@doc """
facies_props_depthbins(fac_ids::AbstractArray{<:Integer,3};
                           bin_m::Real=10.0,
                           n_facies::Int=5,
                           Δz_m::Real)

Compute facies proportions as a function of depth using fixed vertical bins.

# Arguments
- `fac_ids`: 3D array of facies identifiers.
- `bin_m`: Thickness of each depth bin.
- `n_facies`: Number of facies categories.
- `Δz_m`: Vertical grid spacing.

# Returns
- Array of facies proportions per depth bin.

# Algorithm
1. Convert voxel indices to physical depth using `Δz_m`.
2. Partition the vertical domain into bins of size `bin_m`.
3. For each bin:
   - Identify voxels within the bin.
   - Count occurrences of each facies.
4. Normalize counts to proportions.

# Units
- Depth in meters.
- Output is dimensionless proportions.

# Notes
Used for post-processing vertical facies distributions.
""" -> facies_props_depthbins

@doc """
CSV(kwargs...)

Construct a CSV export specification.

# Arguments
- `kwargs...`: Export configuration parameters.

# Returns
- `CSV` specification object.

# Algorithm
1. Wrap keyword arguments in an `IdDict`.
2. Construct `CSV` specification.

# Units
Not applicable.

# Notes
Defines how data are exported to CSV.
\"\"\"
CSV(kwargs...) = CSV(IdDict(kwargs...))



\"\"\"
    header_from_input_exact(input)

Construct a header object directly from model input.

# Arguments
- `input`: Model configuration.

# Returns
- `Header` object.

# Algorithm
1. Compute time axis using `time_axis`.
2. Compute spatial axes using `box_axes`.
3. Construct `Header` with axes and metadata.

# Units
- Distance in meters.
- Time in Myr.

# Notes
Provides exact reconstruction without state dependence.
""" -> header_from_input_exact

@doc """
datacolumn_from_state(state; x_index=1, y_index=1, write_interval=1)

Extract a 1D data column from model state.

# Arguments
- `state`: Model state.
- `x_index`, `y_index`: Spatial indices.
- `write_interval`: Temporal subsampling interval.

# Returns
- `DataColumn` object.

# Algorithm
1. Select column at `(x_index, y_index)`.
2. Subsample time using `write_interval`.
3. Concatenate data along time axis.
4. Construct `DataColumn`.

# Units
- Distance in meters.
- Time in Myr.

# Notes
Used for column-based analysis and export.
""" -> datacolumn_from_state

@doc """
dataslice_from_state_exact(state; y_index=1, write_interval=1)

Extract a 2D slice from model state.

# Arguments
- `state`: Model state.
- `y_index`: Slice index.
- `write_interval`: Temporal subsampling interval.

# Returns
- `DataSlice` object.

# Algorithm
1. Select slice at `y_index`.
2. Subsample time dimension.
3. Concatenate data into slice structure.

# Units
- Distance in meters.
- Time in Myr.

# Notes
Exact reconstruction without interpolation.
""" -> dataslice_from_state_exact

@doc """
datavolume_from_state(state; write_interval=1)

Construct a full 3D data volume from model state.

# Arguments
- `state`: Model state.
- `write_interval`: Temporal subsampling interval.

# Returns
- `DataVolume` object.

# Algorithm
1. Subsample state history.
2. Concatenate fields across time and space.
3. Construct volume representation.

# Units
- Distance in meters.
- Time in Myr.

# Notes
Used for volumetric export.
""" -> datavolume_from_state

@doc """
data_kind(gid::HDF5.Group)

Determine the data type stored in an HDF5 group.

# Arguments
- `gid`: HDF5 group.

# Returns
- Data type identifier.

# Algorithm
1. Inspect group attributes.
2. Parse metadata using `parse_multi_slice`.
3. Return data classification.

# Units
Not applicable.

# Notes
Used for dynamic dispatch during data loading.
""" -> data_kind

@doc """
data_kind(fid::HDF5.File, group)

Determine data type from file and group name.

# Arguments
- `fid`: HDF5 file.
- `group`: Group name.

# Returns
- Data type identifier.

# Algorithm
1. Access group in file.
2. Delegate to `data_kind(gid)`.

# Units
Not applicable.
""" -> data_kind

@doc """
group_datasets(fid::HDF5.File)

Group datasets by type within an HDF5 file.

# Arguments
- `fid`: HDF5 file.

# Returns
- Dictionary mapping dataset types to groups.

# Algorithm
1. Iterate over file keys.
2. Classify each group using `data_kind`.
3. Organize into dictionary.

# Units
Not applicable.
""" -> group_datasets

@doc """
data_header(gid::HDF5.Group)

Reconstruct header metadata from HDF5 group.

# Arguments
- `gid`: HDF5 group.

# Returns
- `DataHeader` object.

# Algorithm
1. Read attributes.
2. Parse axes and metadata.
3. Construct header object.

# Units
Time in Myr.
""" -> data_header

@doc """
read_header(fid)

Read header metadata from an HDF5 file.

# Arguments
- `fid`: HDF5 file.

# Returns
- `Header` object.

# Algorithm
1. Read attributes and axes.
2. Reconstruct header structure.

# Units
- Distance in meters.
- Time in Myr.
""" -> read_header

@doc """
read_data(::Type{Val{dim}}, gid::Union{HDF5.File, HDF5.Group}) where {dim}

Read data of specified dimensionality from HDF5.

# Arguments
- `Val{dim}`: Dimensionality (1, 2, or 3).
- `gid`: File or group.

# Returns
- Data object (`DataColumn`, `DataSlice`, or `DataVolume`).

# Algorithm
1. Parse dataset structure.
2. Extract arrays.
3. Construct appropriate data object.

# Units
- Distance in meters.
- Time in Myr.
""" -> read_data

@doc """
read_data(D::Type{Val{dim}}, filename::AbstractString, group) where {dim}

Read data from file by group name.

# Arguments
- `D`: Dimensionality specifier.
- `filename`: File path.
- `group`: Group name.

# Returns
- Data object.

# Algorithm
1. Open file.
2. Read header.
3. Delegate to `read_data`.

# Units
Inherited from stored data.
""" -> read_data

@doc """
read_volume(args...)

Read 3D volume data.

# Arguments
- `args...`: Passed to `read_data`.

# Returns
- `DataVolume`

# Algorithm
1. Call `read_data(Val{3}, ...)`.
\"\"\"
read_volume(args...) = read_data(Val{3}, args...)



\"\"\"
    read_slice(args...)

Read 2D slice data.

# Arguments
- `args...`: Passed to `read_data`.

# Returns
- `DataSlice`

# Algorithm
1. Call `read_data(Val{2}, ...)`.
\"\"\"
read_slice(args...) = read_data(Val{2}, args...)



\"\"\"
    read_column(args...)

Read 1D column data.

# Arguments
- `args...`: Passed to `read_data`.

# Returns
- `DataColumn`

# Algorithm
1. Call `read_data(Val{1}, ...)`.
\"\"\"
read_column(args...) = read_data(Val{1}, args...)



\"\"\"
    time(header::Header, data::Data)

Return time axis aligned with exported data.

# Arguments
- `header`: Header metadata.
- `data`: Data object.

# Returns
- Time vector.

# Algorithm
1. Subsample `header.axes.t` using `write_interval`.

# Units
Time in Myr.
\"\"\"
time(header::Header, data::Data) = header.axes.t[1:data.write_interval:end]



\"\"\"
    unitful_headers(df::DataFrame)

Return column names including unit annotations.

# Arguments
- `df`: DataFrame.

# Returns
- Vector of strings.

# Algorithm
1. Extract column names and units.
2. Format as strings.

# Units
Not applicable.
\"\"\"
unitful_headers(df::DataFrame) =



\"\"\"
    Unitful.ustrip(df::DataFrame)

Remove units from DataFrame values.

# Arguments
- `df`: DataFrame.

# Returns
- DataFrame with unitless values.

# Algorithm
1. Apply `ustrip` to all columns.
\"\"\"
Unitful.ustrip(df::DataFrame) =



\"\"\"
    write_unitful_csv(io::IO, df::DataFrame)

Write DataFrame to CSV with units in header.

# Arguments
- `io`: Output stream.
- `df`: DataFrame.

# Returns
- `nothing`

# Algorithm
1. Generate unit-aware headers.
2. Strip units from values.
3. Write CSV.
\"\"\"
write_unitful_csv(io, df::DataFrame) =



\"\"\"
    Base.accumulate(f)

Return a curried accumulate function.

# Arguments
- `f`: Reduction function.

# Returns
- Callable function.

# Algorithm
1. Return closure calling `accumulate(f, ...)`.
\"\"\"
Base.accumulate(f) = (args...; kwargs...) -> accumulate(f, args...; kwargs...)



\"\"\"
    age_depth_model(sac::Vector{T}) where {T}

Construct monotonic age-depth model.

# Arguments
- `sac`: Stratigraphic accumulation curve.

# Returns
- Smoothed curve.

# Algorithm
1. Reverse vector.
2. Apply cumulative minimum.
3. Reverse back.

# Units
Distance in meters.
\"\"\"
age_depth_model(sac::Vector{T}) where {T} = sac |> reverse |> accumulate(min) |> reverse



\"\"\"
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
""" -> stratigraphic_column

@doc """
data_export(spec::CSV, header::Header, data)

Export data to CSV.

# Arguments
- `spec`: Export specification.
- `header`: Metadata.
- `data`: Data object.

# Returns
- `nothing`

# Algorithm
1. Iterate over data.
2. Format values.
3. Write to file.
""" -> data_export

@doc """
data_export(::Type{CSVExportTrait{S}}, header::Header, data::DataColumn, label)

Export a single column to CSV.

# Arguments
- Trait type
- `header`
- `data`
- `label`

# Returns
- `nothing`

# Algorithm
1. Validate export configuration.
2. Write column data.
""" -> data_export

@doc """
data_export(::Type{CSVExportTrait{S}}, header::Header, columns)

Export multiple columns.

# Arguments
- Trait type
- `header`
- `columns`

# Returns
- `nothing`

# Algorithm
1. Merge columns.
2. Write combined dataset.
""" -> data_export

@doc """
extract_sac(header::Header, data::DataColumn, label)

Extract stratigraphic accumulation curve.

# Arguments
- `header`
- `data`
- `label`

# Returns
- DataFrame.

# Algorithm
1. Extract thickness vs time.
2. Build DataFrame.
""" -> extract_sac

@doc """
extract_sc(header::Header, data::DataColumn, label)

Extract stratigraphic column.

# Arguments
- `header`
- `data`
- `label`

# Returns
- DataFrame.

# Algorithm
1. Compute stratigraphic column.
2. Convert to table.
""" -> extract_sc

@doc """
extract_wd(header::Header, data::DataColumn, label)

Extract water-depth history.

# Arguments
- `header`
- `data`
- `label`

# Returns
- DataFrame.

# Algorithm
1. Extract water depth.
2. Build table.
""" -> extract_wd

@doc """
initial_state(input::AbstractInput)

Construct the initial model state for the ALCAP model configuration.

# Arguments
- `input`: Model configuration containing grid, forcing, and facies definitions.

# Returns
- State object.

# Algorithm
1. Initialize cellular automaton state using `CellularAutomaton.initial_state`.
2. Perform initial CA update using `CellularAutomaton.step!`.
3. Return initialized state.

# Units
Inherited from input configuration.

# Notes
Extends the base initialization by incorporating CA-based facies initialization.
""" -> initial_state

@doc """
step!(input::Input)

Advance the ALCAP model by one timestep.

# Arguments
- `input`: Model configuration and current state.

# Returns
- `nothing`

# Algorithm
1. Update facies distribution using `CellularAutomaton.step!`.
2. Compute disintegration using `ActiveLayer.disintegrator`.
3. Compute carbonate production using `production`.
4. Transport sediment using `ActiveLayer.transporter`.
5. Update state variables in place.

# Units
- Distance in meters.
- Time in Myr.

# Notes
Defines the full timestep update sequence for the ALCAP model.
""" -> step!

@doc """
write_header(input::AbstractInput, output::AbstractOutput)

Write model-specific metadata for the ALCAP configuration.

# Arguments
- `input`: Model configuration.
- `output`: Output backend.

# Returns
- `nothing`

# Algorithm
1. Iterate over model components.
2. Delegate header writing to component-level implementations using `for_each`.
3. Call `P.write_header` for each component.

# Units
Inherited from component definitions.

# Notes
Extends base header writing with model-specific metadata.
""" -> write_header

@doc """
function initial_state(input::Input)

Implement the `initial_state` step of the current workflow.

# Arguments
- `input::Input`: Model input containing the configuration, forcing functions, grid, time settings, and facies definitions required by this method.

# Returns
- Value constructed by the implementation.

# How it works
1. Allocate the intermediate arrays, plot objects, or output containers needed by the workflow.
2. Delegate lower-level operations to helper routines including `State`.

# Units
This helper relies on the units already carried by its inputs.

# Domain meaning
This method contributes to the core CarboKitten modeling workflow.
""" -> initial_state

@doc """
function step!(input::Input)

Update the target object in place according to the `step!` workflow.

# Arguments
- `input::Input`: Model input containing the configuration, forcing functions, grid, time settings, and facies definitions required by this method.

# Returns
- Nothing. The target state, arrays, or output backend are updated in place.

# How it works
1. Delegate lower-level operations to helper routines including `uniform_production`, `Frame`.

# Units
This helper relies on the units already carried by its inputs.

# Domain meaning
This method contributes to the core CarboKitten modeling workflow.
""" -> step!

@doc """
function write_header(input::AbstractInput, output::AbstractOutput)

Write the `header` product to the configured output backend.

# Arguments
- `input::AbstractInput`: Model input containing the configuration, forcing functions, grid, time settings, and facies definitions required by this method.
- `output::AbstractOutput`: Output backend receiving datasets or attributes.

# Returns
- Value constructed by the implementation.

# How it works
1. Iterate over the relevant grid cells, layers, timesteps, or stored outputs required by the calculation.
2. Delegate lower-level operations to helper routines including `for_each`, `P.write_header`.

# Units
This helper relies on the units already carried by its inputs.

# Domain meaning
This method contributes to the core CarboKitten modeling workflow.
""" -> write_header

@doc """
count_ints(args...)

Count how many arguments are integers.

# Arguments
- `args...`: Arbitrary arguments.

# Returns
- Integer count of arguments of type `Int`.

# Algorithm
1. Recursively traverse arguments.
2. Increment count when encountering `Int`.
3. Return total count.

# Units
Dimensionless.

# Notes
Used to infer dimensionality of slice specifications.
\"\"\"
count_ints(::Int, args...) = 1 + count_ints(args...)
count_ints(_, args...) = count_ints(args...)
count_ints() = 0



\"\"\"
    reduce_slice(s, x, y)

Resolve a slice specification into explicit indices.

# Arguments
- `s`: Slice tuple containing `Int` or `Colon`.
- `x`, `y`: Grid indices.

# Returns
- Tuple of resolved `(x, y)` indices.

# Algorithm
1. Replace `Colon` with provided index.
2. Preserve explicit indices.
3. Return resolved tuple.

# Units
Index space (dimensionless).

# Notes
Used to convert abstract slice definitions into concrete indices.
\"\"\"
reduce_slice(s::Tuple{Colon,Colon}, x, y) = (x, y)
reduce_slice(s::Tuple{Int,Colon}, y::Int) = (s[1], y)
reduce_slice(s::Tuple{Colon,Int}, x::Int) = (x, s[2])



\"\"\"
    parse_slice(s::AbstractString)

Parse a string slice specification.

# Arguments
- `s`: String representing slice (e.g. ":", "10").

# Returns
- Parsed slice element (`Int` or `Colon`).

# Algorithm
1. Interpret ":" as `Colon()`.
2. Otherwise parse integer value.

# Units
Index space.

# Notes
Used for parsing user-defined output specifications.
""" -> parse_slice

@doc """
parse_multi_slice(s::AbstractString)

Parse a multi-dimensional slice specification.

# Arguments
- `s`: Comma-separated slice string.

# Returns
- `Slice2` object.

# Algorithm
1. Split string by comma.
2. Parse each component using `parse_slice`.
3. Wrap result in `Slice2`.

# Units
Index space.
\"\"\"
parse_multi_slice(s::AbstractString) = Slice2(parse_slice.(split(s, ",")))



\"\"\"
    data_kind(args...)

Determine output dimensionality.

# Arguments
- Arguments representing slice specification.

# Returns
- `:column`, `:slice`, or `:volume`.

# Algorithm
1. Count number of integer indices.
2. Classify:
   - 2 integers → `:column`
   - 1 integer → `:slice`
   - 0 integers → `:volume`

# Units
Not applicable.

# Notes
Controls export behavior and data layout.
\"\"\"
data_kind(::Int, ::Int) = :column
data_kind(::Int, _)     = :slice
data_kind(_, ::Int)     = :slice
data_kind(_, _)         = :volume

data_kind(spec::OutputSpec) = data_kind(spec.slice...)



\"\"\"
    new_output

Construct a new output backend.

# Returns
- Output object.

# Algorithm
1. Allocate output structure.
2. Register datasets.
3. Initialize metadata and attributes.

# Units
Depends on stored data.

# Notes
Abstract constructor for output backends.
""" -> new_output

@doc """
add_data_set(out, name, spec)

Register a dataset in the output backend.

# Arguments
- `out`: Output object.
- `name`: Dataset name.
- `spec`: Output specification.

# Returns
- `nothing`

# Notes
Defines what data will be written during simulation.
""" -> add_data_set

@doc """
set_attribute(out, name, value)

Attach metadata attribute to output.

# Arguments
- `out`: Output object.
- `name`: Attribute name.
- `value`: Attribute value.

# Returns
- `nothing`

# Notes
Used for storing metadata (units, axes, parameters).
""" -> set_attribute

@doc """
write_sediment_thickness(out, name, idx, data)

Write sediment thickness to output.

# Arguments
- `out`: Output object.
- `name`: Dataset name.
- `idx`: Time index.
- `data`: Thickness array.

# Returns
- `nothing`

# Algorithm
1. Map timestep index to write index.
2. Store thickness array.

# Units
Distance in meters.

# Notes
Supports 1D, 2D, and 3D outputs.
""" -> write_sediment_thickness

@doc """
write_production(out, name, idx, data)

Write production field to output.

# Arguments
- `out`: Output object.
- `name`: Dataset name.
- `idx`: Time index.
- `data`: Production array.

# Returns
- `nothing`

# Units
m/Myr.

# Notes
First dimension corresponds to facies.
""" -> write_production

@doc """
write_disintegration(out, name, idx, data)

Write disintegration field to output.

# Arguments
- `out`: Output object.
- `name`: Dataset name.
- `idx`: Time index.
- `data`: Disintegration array.

# Returns
- `nothing`

# Units
m/Myr.
""" -> write_disintegration

@doc """
write_deposition(out, name, idx, data)

Write deposition field to output.

# Arguments
- `out`: Output object.
- `name`: Dataset name.
- `idx`: Time index.
- `data`: Deposition array.

# Returns
- `nothing`

# Units
Distance in meters.
""" -> write_deposition

@doc """
state_writer(input, out)

Create a state-level writer.

# Arguments
- `input`: Model configuration.
- `out`: Output backend.

# Returns
- Writer function.

# Algorithm
1. Initialize dataset writers.
2. Return closure for writing state.

# Notes
Handles writing of static or full-state outputs.
""" -> state_writer

@doc """
frame_writer(input, out)

Create a timestep writer function.

# Arguments
- `input`: Model configuration.
- `out`: Output backend.

# Returns
- Function `(idx, state)`.

# Algorithm
1. Return closure that:
   - Checks write interval.
   - Writes all datasets using write_* functions.

# Notes
Core mechanism for writing simulation frames.
""" -> frame_writer

@doc """
facies_percent_from_ids(fac_ids; n_facies)

Compute facies proportions from a categorical voxel block.

# Arguments
- `fac_ids`: 3D array of facies identifiers.
- `n_facies`: Number of facies.

# Returns
- `Dict{Int,Float64}` mapping facies → proportion.

# Algorithm
1. Count occurrences of each facies.
2. Normalize by total voxel count.
3. Return proportions.

# Units
Dimensionless proportions.

# Notes
Used as input for environment classification.
""" -> facies_percent_from_ids

@doc """
classify_block(props, rules; wdepth=NaN, energy=NaN)

Classify a voxel block into an environment category.

# Arguments
- `props`: Facies proportions.
- `rules`: Vector of `EnvironmentRule`.
- `wdepth`: Water depth.
- `energy`: Wave energy.

# Returns
- Environment identifier.

# Algorithm
1. Evaluate each rule using facies proportions and optional physical variables.
2. Apply rule priority order.
3. Return first matching environment.

# Units
- Depth in meters.
- Energy in model units.

# Notes
Core rule-based classification step.
""" -> classify_block

@doc """
env_selector(strength, env_rules, env_names; wdepth=NaN, energy=NaN, nondep_eps=1e-12)

Select environment category from facies strength distribution.

# Arguments
- `strength`: Facies proportions or strengths.
- `env_rules`: Classification rules.
- `env_names`: Ordered environment names.
- `wdepth`: Water depth.
- `energy`: Wave energy.
- `nondep_eps`: Threshold for non-deposition.

# Returns
- Environment index or label.

# Algorithm
1. If total strength < `nondep_eps`, return non-depositional class.
2. Convert inputs to numeric form if unitful.
3. Apply `classify_block`.
4. Map result to `env_names`.

# Units
- Depth in meters.
- Energy in model units.

# Notes
Wrapper around rule-based classification.
""" -> env_selector

@doc """
build_environment_grid(facies_layers, thickness_layers, wdepth_layers, energy_layers;
                           Δi, Δj, Δk, env_rules, n_facies)

Construct a block-wise environment classification grid.

# Arguments
- Layer stacks for facies, thickness, water depth, and energy.
- `Δi`, `Δj`, `Δk`: Block aggregation sizes.
- `env_rules`: Classification rules.
- `n_facies`: Number of facies.

# Returns
- 3D environment grid.

# Algorithm
1. Partition domain into blocks of size `(Δi, Δj, Δk)`.
2. For each block:
   - Compute facies proportions.
   - Compute mean water depth and energy.
   - Classify using `env_selector`.
3. Store result in grid.

# Units
- Distance in meters.
- Output is categorical.

# Notes
Bridges stratigraphic archive → interpreted environment volume.
""" -> build_environment_grid

@doc """
encode_environments(env_grid, env_names)

Convert environment labels to integer codes.

# Arguments
- `env_grid`: Environment grid (labels or symbols).
- `env_names`: Ordered list of environment names.

# Returns
- Integer-encoded grid.

# Algorithm
1. Build mapping `name → index`.
2. Replace values in grid.
3. Return encoded array.

# Units
Categorical indices.
""" -> encode_environments

@doc """
write_environment_block!(fid, encoded, names; Δi, Δj, Δk)

Write environment grid to HDF5.

# Arguments
- `fid`: HDF5 file.
- `encoded`: Environment grid.
- `names`: Environment labels.
- `Δi`, `Δj`, `Δk`: Block sizes.

# Returns
- `nothing`

# Algorithm
1. Create or access HDF5 group.
2. Write encoded grid.
3. Store metadata (names, block sizes).

# Units
Categorical.
""" -> write_environment_block!

@doc """
write_layer_archive!(fid, state)

Write preserved layer archive to HDF5.

# Arguments
- `fid`: HDF5 file.
- `state`: Model state.

# Returns
- `nothing`

# Algorithm
1. Extract layer-wise facies, thickness, and environmental data.
2. Create HDF5 groups.
3. Write datasets.

# Units
- Thickness in meters.
- Facies categorical.

# Notes
Stores full stratigraphic archive.
""" -> write_layer_archive!

@doc """
make_header(input)

Construct HDF5 header metadata.

# Arguments
- `input`: Model configuration.

# Returns
- Header object.

# Algorithm
1. Build axes using `box_axes` and `time_axis`.
2. Include initial topography.
3. Construct header.

# Units
- Distance in meters.
- Time in Myr.
""" -> make_header

@doc """
H5Output(f, input, filename)

Create an HDF5 output backend.

# Arguments
- `f`: Initialization function.
- `input`: Model configuration.
- `filename`: Output file path.

# Returns
- `H5Output` object.

# Algorithm
1. Build header using `make_header`.
2. Open file with `h5open`.
3. Initialize groups.
4. Apply initialization function `f`.

# Units
Inherited from stored data.
""" -> H5Output

@doc """
add_data_set(out::H5Output, name, spec)

Register dataset in HDF5 output.

# Arguments
- `out`: Output backend.
- `name`: Dataset name.
- `spec`: Output specification.

# Returns
- `nothing`

# Algorithm
1. Determine dataset shape using `axis_size`.
2. Create dataset group.
3. Store metadata.
""" -> add_data_set

@doc """
set_attribute(out::H5Output, name, value)

Write attribute to HDF5.

# Arguments
- `out`: Output backend.
- `name`: Attribute path.
- `value`: Attribute value.

# Returns
- `nothing`

# Algorithm
1. Resolve group using `get_group`.
2. Write attribute.
3. Replace existing value if needed.
""" -> set_attribute

@doc """
state_writer(input, out)

Create writer for static state output.

# Returns
- Closure.

# Algorithm
1. Initialize dataset indices.
2. Return function handling state writes.
""" -> state_writer

@doc """
frame_writer(input, out)

Create timestep writer.

# Returns
- Function `(idx, state)`.

# Algorithm
1. Check write interval.
2. Write datasets using write_* methods.
3. Reshape data to match output layout.
""" -> frame_writer

@doc """
run_model(::Type{Model{M}}, input, filename; env_names)

Run model and export HDF5 outputs with environment classification.

# Arguments
- Model type
- `input`: Model configuration
- `filename`: Output file
- `env_names`: Environment labels

# Returns
- Simulation result.

# Algorithm
1. Initialize HDF5 output.
2. Run model.
3. Write layer archive.
4. Build environment grid.
5. Encode environments.
6. Write environment dataset.

# Units
- Distance in meters
- Time in Myr

# Notes
Extends standard run with derived environment export.
""" -> run_model

@doc """
MemoryOutput(input)

Construct an in-memory output backend.

# Arguments
- `input`: Model configuration.

# Returns
- `MemoryOutput` object.

# Algorithm
1. Delegate to `new_output(MemoryOutput, input)`.

# Notes
Stores all output data in memory instead of writing to disk.
\"\"\"
MemoryOutput(input::AbstractInput) = new_output(MemoryOutput, input)



\"\"\"
    new_output(::Type{MemoryOutput}, input)

Initialize memory-based output storage.

# Arguments
- `input`: Model configuration.

# Returns
- `MemoryOutput` object.

# Algorithm
1. Build axes using `time_axis` and `box_axes`.
2. Construct header metadata.
3. Allocate containers for datasets and attributes.

# Units
- Distance in meters
- Time in Myr
""" -> new_output

@doc """
axis_size(sel, n)

Determine size of a selected axis.

# Arguments
- `sel`: Selector (`Colon`, `Int`, or range)
- `n`: Full axis length

# Returns
- Integer size
\"\"\"
axis_size(::Colon, a::Int) = a
axis_size(::Int, _) = 1
axis_size(r::AbstractRange{Int}, _) = length(r)



\"\"\"
    add_data_set(out::MemoryOutput, label, spec)

Register dataset in memory output.

# Arguments
- `out`: Output backend
- `label`: Dataset name
- `spec`: Output specification

# Returns
- `nothing`

# Algorithm
1. Determine dataset dimensionality using `data_kind`.
2. Compute shape using `axis_size`.
3. Allocate array storage.
4. Store metadata (`DataHeader`).
""" -> add_data_set

@doc """
set_attribute(out::MemoryOutput, name, value)

Store metadata attribute in memory.

# Arguments
- `out`: Output backend
- `name`: Attribute name
- `value`: Attribute value

# Returns
- `nothing`

# Notes
Attributes are stored in a dictionary-like structure.
""" -> set_attribute

@doc """
run_model(::Type{Model{M}}, input, output)

Execute a model and write results to an output backend.

# Arguments
- Model type
- `input`: Model configuration
- `output`: Output backend

# Returns
- Final model state

# Algorithm
1. Write header using `M.write_header`.
2. Initialize state using `M.initial_state`.
3. Create writers:
   - `state_writer`
   - `frame_writer`
4. Loop over timesteps:
   - Advance model using `M.step!`
   - Write outputs using frame writer
5. Return final state

# Units
- Distance in meters
- Time in Myr

# Notes
Core execution driver for all output backends.
""" -> run_model

@doc """
patch_stats_3D_from_ids(fac_ids; n_facies=5, r=1)

Compute spatial patch statistics from a 3D facies grid.

# Arguments
- `fac_ids`: 3D array of facies identifiers
- `n_facies`: Number of facies
- `r`: Neighborhood radius

# Returns
- Patch statistics (e.g. connectivity or size metrics)

# Algorithm
1. For each voxel:
   - Inspect neighborhood of radius `r`
2. Identify connected regions of identical facies
3. Compute statistics (e.g. patch size, continuity)
4. Aggregate per facies

# Units
Categorical (facies IDs)

# Notes
Used for post-processing spatial structure of deposits.
""" -> patch_stats_3D_from_ids

@doc """
thickness_from_topk(topk, Δz_m)

Compute platform thickness statistics from voxel top indices.

# Arguments
- `topk`: Top-of-stack voxel index per column.
- `Δz_m`: Vertical voxel size.

# Returns
- `(mean, min, max)` thickness in meters.

# Algorithm
1. Convert voxel indices to thickness: `thickness = topk * Δz_m`.
2. Compute mean, minimum, and maximum across all columns.

# Units
Distance in meters.

# Notes
Provides fast summary metrics of preserved platform thickness.
\"\"\"
thickness_from_topk(topk::AbstractMatrix{<:Integer}, Δz_m::Real) =



\"\"\"
    run_once(input; run_id=1, n_facies=5, nbins=35, return_state=true, env_names)

Run a single simulation and compute summary diagnostics.

# Arguments
- `input`: Model configuration.
- `run_id`: Identifier for the run.
- `n_facies`: Number of facies.
- `nbins`: Number of depth bins.
- `return_state`: Whether to return final state.
- `env_names`: Environment labels.

# Returns
- Summary metrics and optionally final state.

# Algorithm
1. Run model using `run_model`.
2. Extract facies grid and compute proportions.
3. Compute thickness statistics using `thickness_from_topk`.
4. Optionally compute depth-bin summaries.
5. Return results.

# Units
- Distance in meters
- Facies categorical

# Notes
High-level wrapper for single-run analysis.
""" -> run_once

@doc """
run_single_to_csv(input; facies_names, env_names, out_csv, fail_csv, run_id=1, nbins=35, n_facies=5)

Run a simulation and write summary outputs to CSV.

# Arguments
- `input`: Model configuration.
- `facies_names`: Facies labels.
- `env_names`: Environment labels.
- `out_csv`: Output file path.
- `fail_csv`: Failure log path.
- `run_id`: Identifier.
- `nbins`: Number of depth bins.
- `n_facies`: Number of facies.

# Returns
- `nothing`

# Algorithm
1. Execute `run_once`.
2. Flatten results into tabular format.
3. Write successful runs to `out_csv`.
4. Log failures to `fail_csv`.

# Units
- Distance in meters
- Output values dimensionless or categorical

# Notes
Used for batch experiments and parameter sweeps.
""" -> run_single_to_csv

@doc """
flatten_depthbins(props)

Flatten nested depth-bin facies proportions.

# Arguments
- `props`: Nested dictionary `props[bin][facies]`.

# Returns
- Named tuple of scalar values.

# Algorithm
1. Iterate over bins and facies.
2. Generate flattened key-value pairs.
3. Return named tuple.
\"\"\"
flatten_depthbins(props::Dict{Int, Dict{Int, Float64}})



\"\"\"
    flatten_facies(d)

Flatten facies proportions.

# Arguments
- `d`: Dictionary of facies proportions.

# Returns
- Named tuple.

# Algorithm
1. Convert key-value pairs to flat structure.
\"\"\"
flatten_facies(d::Dict{Int, <:Real}) =



# File: src/RunModel.jl

\"\"\"
    run_model(f, ::Type{Model{M}}, input)
    run_model(f, ::Type{Model{M}}, input, state)

Run a model and pass frames to callback `f`.

# Arguments
- `f`: Callback function receiving frames.
- Model type
- `input`: Model configuration
- `state`: Optional pre-initialized state

# Returns
- Final state

# Algorithm
1. Initialize logger.
2. Initialize state if not provided.
3. Loop over timesteps:
   - Call `M.step!`
   - Emit frame to callback `f`
4. Return final state

# Units
Inherited from model.

# Notes
Core execution loop with callback interface.
\"\"\"
run_model(f, ::Type{Model{M}}, input::AbstractInput) where M =



\"\"\"
    run_model(f, ::Type{Model{M}}, input, state)

Run model with explicit initial state.

# Arguments
- `f`: Frame callback
- `input`: Model configuration
- `state`: Initial state

# Returns
- Final state

# Algorithm
1. Use `with_logger` for execution.
2. Iterate for `n_steps(input)`:
   - Update state using `M.step!`
   - Emit frame via `f`
3. Return state
""" -> run_model

@doc """
push_sediment!(col, parcel)

Add sediment parcel to a vertical column.

# Arguments
- `col`: Sediment column (layers × facies).
- `parcel`: Thickness per facies.

# Returns
- `nothing`

# Algorithm
1. Append parcel to top of column.
2. Update layer structure.

# Units
Distance in meters.
""" -> push_sediment!

@doc """
pop_sediment!(col, Δ)

Remove sediment from column.

# Arguments
- `col`: Sediment column.
- `Δ`: Thickness to remove.

# Returns
- Removed sediment vector.

# Algorithm
1. Remove material from top layers.
2. Use `pop_fraction` for partial layers.
""" -> pop_sediment!

@doc """
push_sediment!(sediment, p)

Add sediment to full 3D grid.

# Arguments
- `sediment`: 4D array (facies × x × y × layers).
- `p`: Deposition field.

# Returns
- `nothing`

# Algorithm
1. Loop over grid cells.
2. Add parcel using column-level push.
""" -> push_sediment!

@doc """
peek_sediment(col, Δ)

Inspect removable sediment without modifying column.

# Arguments
- `col`: Sediment column.
- `Δ`: Thickness to inspect.

# Returns
- Sediment vector.

# Algorithm
1. Traverse top layers.
2. Compute removable fraction.
""" -> peek_sediment

@doc """
peek_sediment(sediment, Δ)

Inspect removable sediment for full grid.

# Arguments
- `sediment`: 4D sediment array.
- `Δ`: Thickness.

# Returns
- Sediment field.

# Algorithm
1. Loop over grid.
2. Apply column-level `peek_sediment`.
""" -> peek_sediment

@doc """
pop_sediment!(cols, amount, out)

Remove sediment from grid.

# Arguments
- `cols`: Sediment grid.
- `amount`: Thickness per column.
- `out`: Output array.

# Returns
- `nothing`

# Algorithm
1. Loop over columns.
2. Remove sediment using column-level pop.
3. Store removed material in `out`.
""" -> pop_sediment!

@doc """
pde_stencil(box::Box{BT}, Δt, ν, out, η, C) where {BT<:Boundary{2}}

Assemble the active-layer diffusion stencil for one transport update.

# Arguments
- `box`: Spatial grid and boundary-condition definition.
- `Δt`: Time step used in the stencil update.
- `ν`: Diffusivity or transport coefficient field.
- `out`: Preallocated output array receiving stencil coefficients or update values.
- `η`: Bed elevation or active-layer thickness field used by the stencil.
- `C`: Sediment concentration or transported quantity.

# Returns
- `out`

# Algorithm
1. Traverse the computational stencil implied by `box`.
2. Evaluate the local stencil coefficients from `Δt`, `ν`, geometry, and boundary conditions.
3. Apply `stencil!` to accumulate the local PDE contribution into `out`.
4. Enforce non-negative or bounded values where required.

# Units
- `Δt` in time units used by the transport scheme.
- `ν` in diffusivity units.
- Length-like fields in meters.

# Notes
This routine builds the discrete PDE operator used by active-layer transport.
""" -> pde_stencil

@doc """
advection_coef!(box::Box{BT}, diffusivity, wave_velocity, w, adv, rct) where {BT}

Compute advection and reaction coefficients from a wave-velocity function.

# Arguments
- `box`: Spatial grid and boundary-condition definition.
- `diffusivity`: Diffusivity field or parameter.
- `wave_velocity`: Callable returning local wave-driven transport velocity.
- `w`: Water-depth field.
- `adv`: Preallocated array receiving advection coefficients.
- `rct`: Preallocated array receiving reaction or residual coefficients.

# Returns
- `nothing`

# Algorithm
1. Loop over all grid cells.
2. Evaluate local wave velocity from the water-depth field.
3. Convert velocity and diffusivity into directional advection coefficients.
4. Store the resulting coefficients in `adv` and any residual terms in `rct`.

# Units
- Velocity in the project transport units.
- Distance-like quantities in meters.
- Coefficients in units consistent with the discretized transport equation.

# Notes
This overload is used when wave velocity is supplied as a function of local state.
""" -> advection_coef!

@doc """
advection_coef!(box::Box{BT}, diffusivity,
                    fields::Tuple{AbstractArray,AbstractArray},
                    w, adv, rct) where {BT}

Compute advection and reaction coefficients from precomputed vector fields.

# Arguments
- `box`: Spatial grid and boundary-condition definition.
- `diffusivity`: Diffusivity field or parameter.
- `fields`: Tuple of precomputed transport fields, typically the x- and y-directed velocity components.
- `w`: Water-depth field.
- `adv`: Preallocated array receiving advection coefficients.
- `rct`: Preallocated array receiving reaction or residual coefficients.

# Returns
- `nothing`

# Algorithm
1. Loop over all grid cells.
2. Read the bounded local transport components from `fields`.
3. Combine local transport direction and diffusivity to form advection coefficients.
4. Store the resulting coefficients in `adv` and residual terms in `rct`.

# Units
- Velocity-like fields in the project transport units.
- Distance-like quantities in meters.
- Coefficients in units consistent with the transport discretization.

# Notes
This overload is used when transport fields have already been computed externally.
""" -> advection_coef!

@doc """
max_dt(adv, dx, courant_max)

Compute the maximum stable timestep from the advection field.

# Arguments
- `adv`: Advection coefficient or velocity field.
- `dx`: Horizontal grid spacing.
- `courant_max`: Maximum allowable Courant number for the integration scheme.

# Returns
- Maximum stable timestep.

# Algorithm
1. Compute the largest absolute advection magnitude in `adv`.
2. Apply the Courant stability condition `Δt_max = courant_max * dx / max(|u|)`.
3. Return the resulting timestep limit.

# Units
- `dx` in meters.
- Returned timestep in the time units implied by `adv`.

# Notes
Used to adapt transport timesteps to satisfy numerical stability.
""" -> max_dt

@doc """
transport_dC!(box::Box{BT}, adv, rct, C, dC) where {BT}

Compute the transport tendency for the advected field.

# Arguments
- `box`: Spatial grid and boundary-condition definition.
- `adv`: Advection coefficient field.
- `rct`: Reaction or residual coefficient field.
- `C`: Current transported field.
- `dC`: Preallocated array receiving the tendency.

# Returns
- `nothing`

# Algorithm
1. Loop over all grid cells.
2. Apply the upwind transport update using the local advection coefficients.
3. Add any local reaction or residual contribution from `rct`.
4. Store the net rate of change in `dC`.

# Units
- `C` in concentration or thickness units of the transported quantity.
- `dC` in the corresponding rate units.

# Notes
This routine produces the right-hand side of the transport equation.
""" -> transport_dC!

@doc """
transport!(box::Box{BT}, diffusivity, wave_velocity, C, w, dC) where {BT}

Assemble and evaluate the transport tendency in place.

# Arguments
- `box`: Spatial grid and boundary-condition definition.
- `diffusivity`: Diffusivity field or parameter.
- `wave_velocity`: Wave-driven transport velocity description.
- `C`: Current transported field.
- `w`: Water-depth field.
- `dC`: Preallocated array receiving the transport tendency.

# Returns
- `nothing`

# Algorithm
1. Construct the local transport stencil or coefficients from `diffusivity`, `wave_velocity`, and `w`.
2. Apply the discrete transport operator to `C`.
3. Store the resulting tendency in `dC`.

# Units
- Distance-like quantities in meters.
- Time-like quantities in the transport scheme units.
- `dC` in rate units corresponding to `C`.

# Notes
This is the in-place low-level transport operator used by higher-level wrappers.
""" -> transport!

@doc """
transport(box::Box{BT}, diffusivity, wave_velocity,
              C::AbstractArray{T}, w) where {BT, T}

Allocate and return the transport tendency for a transported field.

# Arguments
- `box`: Spatial grid and boundary-condition definition.
- `diffusivity`: Diffusivity field or parameter.
- `wave_velocity`: Wave-driven transport velocity description.
- `C`: Current transported field.
- `w`: Water-depth field.

# Returns
- Newly allocated array containing the transport tendency.

# Algorithm
1. Allocate an output array with the same shape as `C` and appropriate units.
2. Call `transport!` to compute the transport tendency.
3. Return the populated tendency array.

# Units
- Returned array has rate units consistent with `C` and the transport timestep convention.
- Distance-like inputs are typically in meters.

# Notes
Convenience wrapper around `transport!` for callers that do not provide preallocated storage.
""" -> transport

@doc """
dir(θ)

Convert a propagation angle into a 2D unit direction vector.

# Arguments
- `θ`: Wave propagation angle (radians).

# Returns
- `Vec2(cos(θ), sin(θ))`: Unit vector pointing in the propagation direction.

# Algorithm
1. Evaluate the cosine and sine of the input angle.
2. Construct a 2D vector with x-component `cos(θ)` and y-component `sin(θ)`.

# Units
- `θ` in radians.
- Output is dimensionless.

# Notes
Used to project scalar wave properties into directional transport components.
\"\"\"
dir(θ) = Vec2(cos(θ), sin(θ))



\"\"\"
    _wave_velocity_at_depth(model::WaveModel, pos, t, h, x_fetch)

Compute the wave orbital velocity at a given depth and location.

# Arguments
- `model`: Wave model containing period, amplitude, and attenuation parameters.
- `pos`: Spatial position used to evaluate directional or spatial variation.
- `t`: Simulation time.
- `h`: Local water depth.
- `x_fetch`: Effective fetch controlling wave growth.

# Returns
- 2D velocity vector representing near-bed or depth-dependent wave orbital motion.

# Algorithm
1. Compute wavelength `L` from the wave period using `wavelength_from_period(T, h)`.
2. Compute wave number `k = 2π / L`.
3. Evaluate vertical decay of orbital velocity using an exponential attenuation:
   - Velocity ∝ exp(-k * depth)
4. Scale the velocity magnitude based on fetch-dependent wave growth.
5. Project the scalar velocity onto a direction vector using `dir(...)`.
6. Return the resulting 2D velocity vector.

# Units
- Depth and wavelength in meters.
- Velocity in length per unit time consistent with the model forcing.

# Notes
This function captures the depth decay of wave orbital motion, which controls sediment transport intensity.
""" -> _wave_velocity_at_depth

@doc """
wavelength_from_period(T, h; tol=1e-10, maxiter=100)

Solve the linear wave dispersion relation for wavelength.

# Arguments
- `T`: Wave period.
- `h`: Water depth.
- `tol`: Convergence tolerance for the iterative solver.
- `maxiter`: Maximum number of iterations.

# Returns
- Wavelength `L` satisfying the dispersion relation.

# Algorithm
1. Convert inputs to consistent units if needed (`ustrip`).
2. Solve the dispersion relation:
   - ω² = gk tanh(kh), where ω = 2π / T
3. Use an iterative root-finding scheme:
   a. Initialize wave number `k`.
   b. Update `k` using the dispersion equation.
   c. Evaluate convergence using `abs(Δk) < tol`.
4. Convert wave number to wavelength:
   - `L = 2π / k`
5. Return `L`.

# Units
- `T` in time units.
- `h` and `L` in meters.

# Notes
Handles both shallow- and deep-water limits through the `tanh(kh)` term.
""" -> wavelength_from_period

@doc """
wave_physics_at_cell(model::WaveModel, pos, t, h, x_fetch; dh=0.05u"m")

Evaluate wave-induced hydrodynamic properties at a grid cell.

# Arguments
- `model`: Wave model containing physical parameters.
- `pos`: Spatial position of the cell.
- `t`: Simulation time.
- `h`: Local water depth.
- `x_fetch`: Effective fetch controlling wave development.
- `dh`: Vertical increment used for estimating gradients.

# Returns
- Tuple or structured result containing wave velocity and derived quantities
  such as energy or shear proxies.

# Algorithm
1. Compute orbital velocity at the target depth using `_wave_velocity_at_depth`.
2. Optionally evaluate velocity at a slightly offset depth (`h + dh`) to estimate vertical gradients.
3. Use exponential decay relationships to compute near-bed velocity magnitude.
4. Derive secondary quantities (e.g. energy proxies) from velocity magnitude:
   - Energy ∝ velocity²
5. Apply non-negativity constraints using `max(...)` where needed.
6. Return velocity and derived diagnostics.

# Units
- Depth in meters.
- Velocity in model transport units.
- Derived energy follows project-specific conventions.

# Notes
This routine bridges wave physics and sediment transport by providing physically
consistent forcing fields at each grid cell.
""" -> wave_physics_at_cell

