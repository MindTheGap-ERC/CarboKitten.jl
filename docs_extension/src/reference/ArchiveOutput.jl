"""
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
"""
function write_environment_block! end

"""
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
"""
function write_layer_archive! end

"""
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
"""
function write_sediment_thickness end

"""
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
"""
function write_production end

"""
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
"""
function write_disintegration end

"""
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
"""
function write_deposition end

"""
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
"""
function run_model end

"""
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
"""
function MemoryOutput end