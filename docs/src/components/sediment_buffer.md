# Sediment Buffers

```component-dag
CarboKitten.Components.SedimentBuffer
```
```@docs
CarboKitten.SedimentStack.sediment_layer
```
For our models of transport and denudation it is important to remember the facies of the sediment for some time into the past. One way to do this, is to remember all contributions of sediment in a stack. Every time we transport or erode sediment, we can pop parcels from this stack. In a three-dimensional model, we need a 2d grid of stacks. Each stack would have its own memory management (which is computationally expensive), and most resources are spent on areas with very little accretion. (In fact, this is what CarboCAT does. We believe this design choice is the main contributor to the difference in run-time between CarboCAT and CarboKitten).

Instead, we choose a fixed size sediment buffer. Each cell in the buffer represents a parcel of sediment, where we store the relative fractions of each contributing facies. This buffer is only used to determine the facies of disintegrated sediment. The output of the overal model is still the amount of sediment for each iteration.

We can't stress this enough: any inaccuracy in using a fixed size buffer with a chosen granularity only impacts the precision of the composition of transported sediment. Even then, the schema is conservative: no sediment is lost unless erosion is so rampant that it eats through the entire sediment stack. In that case, a simulation should be run with a larger buffer.

## Data structure

While the sediment buffer is allocated as a single 4-dimensional array (depth, facies, $x$, $y$), it is best to explain its functioning from the perspective of a single cell in our model. We are left with two dimensions: depth (rows) and facies (columns).

We choose to have the head of our sediment stack always be at the first row. When sediment out-grows the buffer, the deepenst layers are dropped from memory. The head can contain an incomplete amount of sediment, while all rows below the head are either full or empty. When sediment is pushed to the stack and the head row overflows, all rows are copied down one row and the surplus is assigned to the now empty head row. The inverse happens when removing (popping) material from the stack (in computer science stacks are pushed on and popped from). This process is illustrated below.

![sediment buffer diagram](../fig/sediment-buffer.svg)

Above we see a buffer. First we push a parcel of size $3/4$, then we pop an amount of $1/2$. This popped parcel will have different fractions from the pushed one, since it also draws from the half filled row that was in the stack before pushing. In this sense, a small amount of facies mixing will take place, depending on the depositional resolution chosen.

Our implementation is such that each cell in the buffer is contiguous in memory. Thus, copying rows of unstrided memory should be very efficient, although the performance remains to be tested.

## Implementation

We define two functions `push_sediment!` and `pop_sediment!`. Given a $s \times n$ matrix, where $n$ is the number of facies types and $s$ is the depth of the stack, we can grow and shrink sediment. These functions are unit-free, setting $\Delta z$ to be equal to 1.

``` {.julia file=test/SedimentStackSpec.jl}
@testset "SedimentStack" begin
  using CarboKitten.SedimentStack: push_sediment!, pop_sediment!
  stack = zeros(Float64, 10, 3)
  @test pop_sediment!(stack, 0.0) == [0.0, 0.0, 0.0]
  push_sediment!(stack, [5.0, 0, 0])
  @test pop_sediment!(stack, 1.5) == [1.5, 0.0, 0.0]
  push_sediment!(stack, [0.0, 2.0, 0.0])   # (0 0.5) (0 1) (0.5 0.5) (1 0) ...
  @test pop_sediment!(stack, 2.0) == [0.25, 1.75, 0.0]
  @test pop_sediment!(stack, 1.5) == [1.25, 0.25, 0.0]
  @test pop_sediment!(stack, 0.0) == [0.0, 0.0, 0.0]
end

@testset "SedimentArray" begin
  using CarboKitten.SedimentStack: push_sediment!, peek_sediment
  sediment = zeros(Float64, 10, 3, 5, 5)
  for x in 1:10
    production = rand(3, 5, 5)
    push_sediment!(sediment, production)
  end
  a = peek_sediment(sediment, 1.0)
  @test all(sum(a; dims=1) .≈ 1.0)
end
@testset "Sediment layer reconstruction" begin
    using CarboKitten.SedimentStack: sediment_layer

    deposition = zeros(Float64, 3, 1, 1, 4)
    disintegration = zeros(Float64, 3, 1, 1, 4)

    deposition[1, 1, 1, 1] = 1.0
    deposition[2, 1, 1, 2] = 1.0
    deposition[3, 1, 1, 3] = 1.0

    top, present = sediment_layer(deposition, disintegration, 3)
    @test present[1, 1]
    @test top[:, 1, 1] ≈ [0.0, 0.0, 1.0]

    middle, present = sediment_layer(deposition, disintegration, 3; depth=1.0)
    @test present[1, 1]
    @test middle[:, 1, 1] ≈ [0.0, 1.0, 0.0]

    bottom, present = sediment_layer(deposition, disintegration, 3; depth=2.0)
    @test present[1, 1]
    @test bottom[:, 1, 1] ≈ [1.0, 0.0, 0.0]

    combined, present = sediment_layer(
        deposition,
        disintegration,
        3;
        thickness=2.0,
    )
    @test present[1, 1]
    @test combined[:, 1, 1] ≈ [0.0, 1.0, 1.0]

    earlier, present = sediment_layer(deposition, disintegration, 2)
    @test present[1, 1]
    @test earlier[:, 1, 1] ≈ [0.0, 1.0, 0.0]

    absent, present = sediment_layer(deposition, disintegration, 3; depth=3.0)
    @test !present[1, 1]
    @test iszero(sum(absent[:, 1, 1]))

    eroded_deposition = zeros(Float64, 2, 1, 1, 3)
    eroded_disintegration = zeros(Float64, 2, 1, 1, 3)
    eroded_deposition[1, 1, 1, 1] = 1.0
    eroded_deposition[2, 1, 1, 2] = 1.0
    eroded_disintegration[2, 1, 1, 3] = 1.0

    exposed, present = sediment_layer(
        eroded_deposition,
        eroded_disintegration,
        3,
    )
    @test present[1, 1]
    @test exposed[:, 1, 1] ≈ [1.0, 0.0]

    # A finite layer thickness combines a very small final event with the
    # underlying sediment instead of letting that event dominate the map.
    thin_deposition = zeros(Float64, 2, 1, 1, 2)
    thin_disintegration = zeros(Float64, 2, 1, 1, 2)
    thin_deposition[1, 1, 1, 1] = 1.0
    thin_deposition[2, 1, 1, 2] = 0.01

    smoothed, present = sediment_layer(
        thin_deposition,
        thin_disintegration,
        2;
        thickness=1.0,
    )
    @test present[1, 1]
    @test smoothed[:, 1, 1] ≈ [0.99, 0.01]

    # Coverage must be nested with depth: once a preserved column is too thin
    # at one depth, it must remain absent at every greater depth.
    shallow, present_shallow = sediment_layer(
        deposition,
        disintegration,
        3;
        depth=3.0,
    )
    deep, present_deep = sediment_layer(
        deposition,
        disintegration,
        3;
        depth=5.0,
    )
    @test !present_shallow[1, 1]
    @test !present_deep[1, 1]
    @test iszero(sum(shallow))
    @test iszero(sum(deep))
end
```

### Pushing sediment

The single-cell version of `push_sediment!` takes as argument `col` a column (physically speaking a column of sediment) represented by a $s \times n$-matrix and a parcel a $n$-vector.

``` {.julia #sediment-stack-impl}
function push_sediment!(col::AbstractMatrix{F}, parcel::AbstractVector{F}) where F <: Real
    @assert size(col, 2) == length(parcel) "column $(size(col)) doesn't match parcel $(size(parcel))"
    <<push-sediment>>
end
```

First we check if the amount of sediment is larger than the buffer. A warning is printed and the entire buffer filled in the specified fractions.

``` {.julia #push-sediment}
mass = sum(parcel)
if mass == 0.0
    return
end

if mass > size(col)[1]
    @warn "pushing a very large parcel of sediment: $mass times depositional resolution"
    frac = parcel ./ mass
    col .= frac
    return
end
```

First we determine the total sediment amount $\Delta$, being the sum of the parcel, as well as the amount of sediment in our *bucket*, the head row.

``` {.julia #push-sediment}
bucket = sum(col[1, :])
@assert bucket >= 0.0 && bucket <= 1.0
```

If the bucket has enough space left for the parcel, we can just add the parcel to the bucket and return.

``` {.julia #push-sediment}
if bucket + mass < 1.0
    col[1,:] .+= parcel
    return
end
```

Otherwise, we compute the normalized fractions `frac` of facies in the parcel. We add as much sediment as we can to fill the bucket and copy rows down as far as needed.

``` {.julia #push-sediment}
frac = parcel ./ mass
col[1,:] .+= frac .* (1.0 - bucket)
mass -= (1.0 - bucket)
n = floor(Int64, mass)

col[n+2:end,:] .= col[1:end-n-1,:]
```

If the parcel has enough material left to fill more rows, those are all filled with the fractions in `frac`. The head row is assigned whatever is left.

``` {.julia #push-sediment}
na = [CartesianIndex()]
col[2:n+1,:] .= frac[na,:]
mass -= n
col[1,:] .= frac .* mass
```

### Popping sediment

Similar to `push_sediment!` we have `pop_sediment!`. We give `pop_sediment!` the sedimentary column `col` and the total amount of sediment we require. There is a bit that we will reuse called `pop_fraction`, which only works if the amount of popped sediment is lower than the contents of the bucket.

``` {.julia #sediment-stack-impl}
@inline function pop_fraction(col::AbstractMatrix{F}, mass::F) where F <: Real
    bucket = sum(col[1,:])
    if mass == 0 || bucket == 0
        return zeros(F, size(col)[2])
    end

    @assert mass < bucket "pop_fraction can only pop from the top cell: $(col), $(mass)"
    parcel = (mass / bucket) .* col[1,:]
    col[1,:] .-= parcel
    return parcel
end

function pop_sediment!(col::AbstractMatrix{F}, Δ::F) where F <: Real  # -> Vector{F}
    <<pop-sediment>>
end
```

We start by computing the bucket size again. If it is greater than the required amount, we call `pop_fraction`.

``` {.julia #pop-sediment}
bucket = sum(col[1,:])
@assert bucket >= 0.0

if Δ < bucket
  return pop_fraction(col, Δ)
end
```

Otherwise, we start a parcel with the contents of the bucket. Add to that the remaining material in rows below. Now we copy rows from below, setting the bottom $n$ rows to 0. The last step is to call `pop_fraction` one more time with the remaining required amount.

``` {.julia #pop-sediment}
parcel = copy(col[1,:])
Δ -= bucket
n = floor(Int64, Δ)

if n > (size(col)[1] - 2)
    @error "too much material popped of the stack: Δ = $Δ"
    parcel .+= sum(col; dims=1)'
    col .= 0.0
    return parcel
end

parcel .+= sum(col[2:n+1,:]; dims=1)'
col[1:end-n-1, :] = col[n+2:end, :]
col[end-n-1:end, :] .= 0
Δ -= n

parcel .+= pop_fraction(col, Δ)
return parcel
```

### Peeking

Instead of popping sediment, we can also *peek* at the stack with `peek_sediment!`, which is a non-destructive way to inspect what the returned parcel would be if we were to call `pop_sediment!` with the same arguments.

```@raw html
<details><summary>SedimentStack impl</summary>
```

``` {.julia file=src/SedimentStack.jl}
module SedimentStack

export push_sediment!, pop_sediment!, peek_sediment, sediment_layer
<<sediment-stack-impl>>

function push_sediment!(sediment::AbstractArray{F, 4}, p::AbstractArray{F, 3}) where F <: Real
  _, x, y = size(p)
  @assert size(sediment, 2) == size(p, 1) "shapes for sediment $(size(sediment)) doesn't match p $(size(p))"
  @views for i in CartesianIndices((x, y))
    push_sediment!(sediment[:, :, i[1], i[2]], p[:, i[1], i[2]])
  end
end

function peek_sediment(col::AbstractMatrix{F}, Δ::F) where F <: Real  # -> Vector{F}
  if Δ == 0
      return zeros(F, size(col)[2])
  end

  bucket = sum(col[1,:])
  if Δ < bucket
    parcel = (Δ / bucket) .* col[1,:]
    return parcel
  end

  parcel = copy(col[1,:])
  Δ -= bucket
  n = floor(Int64, Δ)

  parcel .+= sum(col[2:n+1,:]; dims=1)'
  Δ -= n

  last_bit = (Δ / sum(col[n+2,:])) .* col[n+2,:]
  parcel .+= last_bit

  return parcel
end

function peek_sediment(sediment::AbstractArray{F,4}, Δ::F) where F <: Real
  _, f, x, y = size(sediment)
  out = Array{F, 3}(undef, f, x, y)
  for i in CartesianIndices((x, y))
    out[:, i[1], i[2]] = peek_sediment(@view(sediment[:, :, i[1], i[2]]), Δ)
  end
  return out
end

"""
    sediment_layer(deposition, disintegration, time_index;
                   depth=0.0, thickness=1.0,
                   amount_to_cells=Float64)

Reconstruct the sediment stack through `time_index` and return a preserved
interval below the sediment surface.

`deposition` and `disintegration` have dimensions `(facies, x, y, time)`.
`depth` and `thickness` are expressed in sediment-buffer cells. `depth` is the
thickness removed below the local top of every reconstructed preserved column.
`amount_to_cells` converts one sediment amount to the same dimensionless units.
The result is `(layer, present)`, where `layer` has dimensions
`(facies, x, y)` and `present` marks cells that contain the requested interval.
"""
function sediment_layer(
    deposition::AbstractArray{T,4},
    disintegration::AbstractArray{T,4},
    time_index::Integer;
    depth::Real = 0.0,
    thickness::Real = 1.0,
    amount_to_cells = Float64,
) where T
    size(deposition) == size(disintegration) ||
        throw(DimensionMismatch(
            "deposition and disintegration must have the same shape",
        ))

    n_facies, nx, ny, n_times = size(deposition)
    1 <= time_index <= n_times ||
        throw(ArgumentError("time_index must be between 1 and $(n_times)"))
    depth >= 0.0 || throw(ArgumentError("depth must be non-negative"))
    thickness > 0.0 || throw(ArgumentError("thickness must be positive"))

    function parcel_mass(data, i, j, k)
        mass = 0.0
        @inbounds for f in 1:n_facies
            amount = Float64(amount_to_cells(data[f, i, j, k]))
            amount >= 0.0 ||
                throw(ArgumentError("sediment amounts must be non-negative"))
            mass += amount
        end
        return mass
    end

    # Determine only the stack depth needed to preserve the requested final
    # interval. In the absence of erosion this is approximately
    # depth + thickness, rather than the complete accumulated succession.
    required_capacity = depth + thickness
    @inbounds for j in 1:ny, i in 1:nx
        height = 0.0
        maximum_height = 0.0

        for k in 1:time_index
            height = max(
                0.0,
                height - parcel_mass(disintegration, i, j, k),
            )
            height += parcel_mass(deposition, i, j, k)
            maximum_height = max(maximum_height, height)
        end

        if height > depth
            sample_bottom = max(0.0, height - depth - thickness)
            required_capacity = max(
                required_capacity,
                maximum_height - sample_bottom,
            )
        end
    end

    # Extra empty rows keep complete-buffer pops away from the edge case in
    # pop_sediment! while leaving the reconstructed result unchanged.
    n_layers = max(3, ceil(Int, required_capacity) + 3)

    layer = zeros(Float64, n_facies, nx, ny)
    present = falses(nx, ny)
    column = zeros(Float64, n_layers, n_facies)
    parcel = zeros(Float64, n_facies)

    @inbounds for j in 1:ny, i in 1:nx
        fill!(column, 0.0)
        available = 0.0

        for k in 1:time_index
            eroded = min(
                parcel_mass(disintegration, i, j, k),
                available,
            )

            if eroded > 0.0
                pop_sediment!(column, eroded)
                available -= eroded
            end

            deposited = 0.0
            for f in 1:n_facies
                amount = Float64(amount_to_cells(deposition[f, i, j, k]))
                amount >= 0.0 ||
                    throw(ArgumentError("sediment amounts must be non-negative"))
                parcel[f] = amount
                deposited += amount
            end

            if deposited > 0.0
                push_sediment!(column, parcel)
                available = min(Float64(n_layers), available + deposited)
            end
        end

        available <= depth && continue

        if depth > 0.0
            pop_sediment!(column, depth)
            available -= depth
        end

        sampled = min(thickness, available)
        sampled <= 0.0 && continue

        layer[:, i, j] .= pop_sediment!(column, sampled)
        present[i, j] = true
    end

    return layer, present
end

function pop_sediment!(cols::AbstractArray{F, 4}, amount::AbstractArray{F, 2}, out::AbstractArray{F, 3}) where F <: Real
  @views for i in CartesianIndices(amount)
      out[:, i[1], i[2]] = pop_sediment!(cols[:, :, i[1], i[2]], amount[i[1], i[2]])
  end
end

end # module
```

```@raw html
</details>
```

### Extracting a stratigraphic layer

Map views reconstruct stratigraphy by replaying deposition and disintegration
through the same sediment-stack operations used by the model. After the stack
has been built to the requested time, `sediment_layer` removes the selected
depth below the local top of each preserved column and returns a finite
interval. The finite thickness prevents a very small last sedimentation event
from producing unstable facies fractions.

The helper is unit-free, like the rest of `SedimentStack`. A caller working in
physical lengths converts sediment amounts, depth, and layer thickness with the
model's depositional resolution.

## Component

``` {.julia file=src/Components/SedimentBuffer.jl}
@compose module SedimentBuffer
@mixin Boxes, FaciesBase, WaterDepth

using StaticArrays
using Unitful

using ..Common
using CarboKitten.SedimentStack: pop_sediment!, push_sediment!, peek_sediment

export pop_sediment, push_sediment

@kwdef struct Input <: AbstractInput
    sediment_buffer_size::Int = 50
    depositional_resolution::Amount = 0.5u"m"
end

@kwdef mutable struct State <: AbstractState
    sediment_buffer::Array{Float64,4}
    sediment_thickness::Array{Height, 2}
end

@constructor _initial_state(input)::State[sediment_buffer, sediment_thickness] = (
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...),
    sediment_thickness = zeros(Amount, input.box.grid_size...))

function push_sediment(input::AbstractInput)
    res = input.depositional_resolution
    n_f = n_facies(input)
    n_g = input.box.grid_size

    function (state::AbstractState, sediment::Array{Amount, 3})
        for i in CartesianIndices(n_g)
            total = sum(@view sediment[:, i[1], i[2]])
            state.sediment_thickness[i] += total
            state.bathymetry[i] += total
            v = SVector{n_f, Float64}(sediment[:, i[1], i[2]] ./ res .|> NoUnits)
            push_sediment!(view(state.sediment_buffer, :, :, i[1], i[2]), v)
        end
    end
end

function pop_sediment(input::AbstractInput)
    res = input.depositional_resolution
    n_g = input.box.grid_size

    function (state::AbstractState, amount::Array{Amount, 2}, out::Array{Amount, 3})
        for i in CartesianIndices(n_g)
            state.sediment_thickness[i] -= amount[i]
            state.bathymetry[i] -= amount[i]
            v = pop_sediment!(@view(state.sediment_buffer[:, :, i[1], i[2]]), amount[i] ./ res .|> NoUnits)
            view(out, :, i[1], i[2]) .= v .* res
        end
    end
end

end
```

### Component Test

``` {.julia file=test/Components/SedimentBufferSpec.jl}
using CarboKitten
using CarboKitten.Components.Common: Amount
import CarboKitten.Components.SedimentBuffer as SB

@testset "Components/SedimentBuffer" begin
    input = SB.Input(
        box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
        time = TimeProperties(Δt=1.0u"yr", steps=10),
        facies = [SB.Facies()],
        sediment_buffer_size = 10,
        depositional_resolution = 1.0u"m")
    state = SB._initial_state(input)

    @test size(state.sediment_thickness) == (10, 1)
    @test size(state.sediment_buffer) == (10, 1, 10, 1)

    push! = SB.push_sediment(input)
    pop! = SB.pop_sediment(input)

    push!(state, reshape((1:10) .* 0.5u"m" |> collect, (1, 10, 1)))
    @test state.sediment_buffer[1, 1, :, 1] == repeat([0.5, 0.0], 5)
    # no initial topography, no subsidence
    @test state.sediment_thickness ≈ state.bathymetry

    buffer = zeros(Amount, 1, 10, 1)
    pop!(state, state.sediment_thickness ./ 2, buffer)
    @test reshape(state.sediment_thickness, (1, 10, 1)) ≈ buffer
    @test state.sediment_buffer[1, 1, :, 1] .% 1.0 ≈ repeat([0.25, 0.5, 0.75, 0.0], 3)[1:10]
    # no initial topography, no subsidence
    @test state.sediment_thickness ≈ state.bathymetry
end
```
