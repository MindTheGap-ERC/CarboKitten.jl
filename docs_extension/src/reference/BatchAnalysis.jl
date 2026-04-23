"""
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
"""
function patch_stats_3D_from_ids end

"""
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
"""
thickness_from_topk(topk::AbstractMatrix{<:Integer}, Δz_m::Real) = nothing



"""
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
"""
function run_once end

"""
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
"""
function run_single_to_csv end

"""
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
"""
flatten_depthbins(props::Dict{Int, Dict{Int, Float64}}) = nothing



"""
    flatten_facies(d)

Flatten facies proportions.

# Arguments
- `d`: Dictionary of facies proportions.

# Returns
- Named tuple.

# Algorithm
1. Convert key-value pairs to flat structure.
"""
flatten_facies(d::Dict{Int, <:Real}) = nothing