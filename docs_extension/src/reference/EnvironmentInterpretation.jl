"""
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
- `(props, edges)` where:
  - `props` is a matrix of facies proportions per depth bin.
  - `edges` is the vector of depth-bin edges.

# Algorithm
1. Convert voxel indices to physical depth using `Δz_m`.
2. Partition the vertical domain into bins of size `bin_m`.
3. For each bin:
   - Identify voxels within the bin.
   - Count occurrences of each facies.
4. Normalize counts to proportions within each bin.

# Units
- Depth in meters.
- Output is dimensionless proportions.

# Notes
This function is occurrence-based: proportions are computed from facies voxel counts
within each depth bin, not from depositional thickness.
For thickness-weighted depth-bin summaries, use `facies_props_depthbins_thickness`.
"""
function facies_props_depthbins end

"""
    facies_props_depthbins_thickness(
        facies_layers::AbstractVector{<:AbstractMatrix{<:Integer}},
        thickness_layers::AbstractVector{<:AbstractMatrix{<:Real}};
        bin_m::Real=10.0,
        n_facies::Int=5,
        Δz_m::Real
    )

Compute facies proportions as a function of depth using fixed vertical bins,
weighted by thickness contribution from archived layer stacks.

# Arguments
- `facies_layers`: Vector of 2D facies matrices, one per deposited layer.
- `thickness_layers`: Vector of 2D thickness matrices with the same structure as `facies_layers`.
- `bin_m`: Thickness of each depth bin.
- `n_facies`: Number of facies categories.
- `Δz_m`: Vertical grid spacing used to map layer index to depth.

# Returns
- `(props, edges)` where:
  - `props` is a `nbins × n_facies` matrix of facies percentages by thickness.
  - `edges` is the vector of depth-bin edges in meters.

# Algorithm
1. Convert layer index to physical depth using `Δz_m`.
2. Partition the vertical domain into bins of size `bin_m`.
3. For each depth bin:
   - Identify the layer indices whose centers fall inside the bin.
   - Sum thickness contribution for each facies across those layers.
4. Normalize facies thickness by total thickness in the bin.
5. Convert to percentages.

# Units
- Depth in meters.
- Input thickness in model thickness units, typically meters.
- Output in percent of total thickness per depth bin.

# Notes
This is the thickness-weighted analogue of `facies_props_depthbins`.
Unlike the voxel-count version, this function uses archived deposited thickness
and therefore preserves relative facies contribution by thickness.
"""
function facies_props_depthbins_thickness end

"""
    facies_thickness_from_layers(facies_layers, thickness_layers; n_facies=5)

Accumulate thickness contribution by facies from archived layer stacks.

# Arguments
- `facies_layers`: Stack of facies identifiers for deposited layers.
- `thickness_layers`: Stack of layer thicknesses with the same shape as `facies_layers`.
- `n_facies`: Number of facies categories.

# Returns
- `Dict{Int,Float64}` mapping facies → total thickness contribution.

# Algorithm
1. Iterate through facies and thickness layers together.
2. For each deposited layer element:
   - Read its facies identifier.
   - Add its thickness to the corresponding facies total.
3. Return accumulated thickness by facies.

# Units
- Input thickness in model thickness units, typically meters.
- Output in the same thickness units.

# Notes
Unlike voxel-count summaries, this function is thickness-weighted and preserves
the volumetric contribution of each facies in the stratigraphic archive.
"""
function facies_thickness_from_layers end

"""
    facies_percent_from_ids(fac_ids; n_facies=5)

Compute facies proportions from a categorical voxel block.

# Arguments
- `fac_ids`: 3D array of facies identifiers.
- `n_facies`: Number of facies categories.

# Returns
- `Dict{Int,Float64}` mapping facies → percentage contribution by voxel count.

# Algorithm
1. Count occurrences of each facies ID in the voxel block.
2. Sum the total number of occupied voxels.
3. Normalize facies counts by the total.
4. Convert to percentages.

# Units
- Output in percent of total occupied voxels.

# Notes
This function is occurrence-based, not thickness-based.
Use `facies_percent_from_layers` when proportions should represent relative
thickness contribution from archived layer stacks.
"""
function facies_percent_from_ids end

"""
    facies_percent_from_layers(facies_layers, thickness_layers; n_facies=5)

Compute facies proportions from archived layer stacks using thickness weighting.

# Arguments
- `facies_layers`: Stack of facies identifiers for deposited layers.
- `thickness_layers`: Stack of layer thicknesses with the same shape as `facies_layers`.
- `n_facies`: Number of facies categories.

# Returns
- `Dict{Int,Float64}` mapping facies → percentage contribution by thickness.

# Algorithm
1. Accumulate total thickness contribution for each facies.
2. Sum thickness over all facies.
3. Normalize each facies thickness by the total thickness.
4. Convert to percentages.

# Units
- Output in percent of total deposited thickness.

# Notes
This function is thickness-based, not occurrence-based.
It should be used when facies proportions are intended to represent relative
thickness contribution rather than voxel abundance.
"""
function facies_percent_from_layers end

"""
    classify_block(props, rules; wdepth=NaN, energy=NaN)

Classify a block into an environment category from facies proportions and optional
physical covariates.

# Arguments
- `props`: Facies proportions or relative facies strengths for the block.
- `rules`: Vector of `EnvironmentRule`.
- `wdepth`: Water depth.
- `energy`: Wave energy.

# Returns
- Environment identifier.

# Algorithm
1. Evaluate each rule against the supplied facies proportions and optional physical variables.
2. Apply rules in priority order.
3. Return the first matching environment.

# Units
- Facies input is dimensionless.
- Depth in meters if supplied.
- Energy in model units if supplied.

# Notes
`props` may come from either occurrence-based or thickness-based facies summaries,
depending on the upstream workflow. The classifier itself only consumes the relative
facies distribution it is given.
"""
function classify_block end

"""
    env_selector(strength, env_rules, env_names; wdepth=NaN, energy=NaN, nondep_eps=1e-12)

Select an environment category from a facies distribution and optional physical variables.

# Arguments
- `strength`: Facies proportions or facies strengths.
- `env_rules`: Classification rules.
- `env_names`: Ordered environment names.
- `wdepth`: Water depth.
- `energy`: Wave energy.
- `nondep_eps`: Threshold for non-deposition.

# Returns
- Environment index or label.

# Algorithm
1. If total facies strength is below `nondep_eps`, assign the non-depositional class.
2. Convert unitful values to plain numeric values if needed.
3. Apply `classify_block`.
4. Map the result to `env_names`.

# Units
- Facies input is dimensionless.
- Depth in meters if supplied.
- Energy in model units if supplied.

# Notes
This is a wrapper around rule-based classification and is agnostic to whether the
input facies distribution was computed from voxel occurrence or thickness contribution.
"""
function env_selector end

"""
    build_environment_grid(facies_layers, thickness_layers, wdepth_layers, energy_layers;
                           Δi, Δj, Δk, env_rules, n_facies)

Construct a block-wise environment classification grid from archived stratigraphic layers.

# Arguments
- `facies_layers`: Layered facies archive.
- `thickness_layers`: Layered thickness archive.
- `wdepth_layers`: Layered water-depth archive.
- `energy_layers`: Layered energy archive.
- `Δi`, `Δj`, `Δk`: Block aggregation sizes.
- `env_rules`: Classification rules.
- `n_facies`: Number of facies.

# Returns
- 3D environment grid.

# Algorithm
1. Partition the domain into blocks of size `(Δi, Δj, Δk)`.
2. For each block:
   - Compute facies proportions from the archived facies/thickness layers.
   - Compute representative water depth and energy.
   - Classify the block using `env_selector`.
3. Store the resulting environment class in the output grid.

# Units
- Thickness and depth in meters if the archive is in physical units.
- Output is categorical.

# Notes
When paired with `facies_percent_from_layers`, the environment grid is based on
thickness-weighted facies composition rather than simple voxel occurrence.
This bridges the stratigraphic archive to an interpreted environment volume.
"""
function build_environment_grid end

"""
    encode_environments(env_grid, env_names)

Convert environment labels to integer codes.

# Arguments
- `env_grid`: Environment grid containing labels, symbols, or names.
- `env_names`: Ordered list of environment names.

# Returns
- Integer-encoded grid.

# Algorithm
1. Build the mapping `name → index`.
2. Replace each environment value in `env_grid` with its corresponding integer code.
3. Return the encoded array.

# Units
Categorical indices.

# Notes
This is a post-processing utility and does not depend on how facies proportions
were originally computed.
"""
function encode_environments end