
# Reconstruction and state-to-data conversion
"""
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
""" 
function build_vertical_coordinates end 

"""
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
"""
function build_selected_surfaces end

"""
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
"""
function resample_to_regular_grid end

"""
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
"""
function thickness_field end 

"""
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
"""
function aligned_time end

"""
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
"""
function elevation end

"""
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
"""
function elevation_physical end

"""
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
"""
function header_from_input_exact end

"""
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
"""
function datacolumn_from_state end

"""
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
"""
function dataslice_from_state_exact end

"""
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
"""
function datavolume_from_state end