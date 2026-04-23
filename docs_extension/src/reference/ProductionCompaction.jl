"""
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
"""
function production end


"""
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
"""
function CompactionState end


"""
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
"""
function depth_multiplier end

"""
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
"""
function time_multiplier end

"""
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
"""
function porosity_at_depth end

"""
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
"""
function compacted_thickness end

"""
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
"""
function compact_column! end

"""
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
"""
function update_sediment_height_from_compaction! end

"""
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
"""
function merge_layers! end

"""
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
"""
function extract_compacted_layer_maps! end

"""
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
"""
function project_compacted_layers_to_cube! end

"""
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
"""
function bury_deposition! end

"""
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
"""
production_rate(f, w, t) = nothing



"""
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
"""
function capped_production end

"""
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
"""
function uniform_production end

"""
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
"""
function initial_topography end

"""
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
"""
function water_depth end