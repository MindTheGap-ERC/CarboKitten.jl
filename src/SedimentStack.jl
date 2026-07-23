# ~/~ begin <<docs/src/components/sediment_buffer.md#src/SedimentStack.jl>>[init]
module SedimentStack

export push_sediment!, pop_sediment!, peek_sediment, sediment_layer
# ~/~ begin <<docs/src/components/sediment_buffer.md#sediment-stack-impl>>[init]
function push_sediment!(col::AbstractMatrix{F}, parcel::AbstractVector{F}) where F <: Real
    @assert size(col, 2) == length(parcel) "column $(size(col)) doesn't match parcel $(size(parcel))"
    # ~/~ begin <<docs/src/components/sediment_buffer.md#push-sediment>>[init]
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
    # ~/~ end
    # ~/~ begin <<docs/src/components/sediment_buffer.md#push-sediment>>[1]
    bucket = sum(col[1, :])
    @assert bucket >= 0.0 && bucket <= 1.0
    # ~/~ end
    # ~/~ begin <<docs/src/components/sediment_buffer.md#push-sediment>>[2]
    if bucket + mass < 1.0
        col[1,:] .+= parcel
        return
    end
    # ~/~ end
    # ~/~ begin <<docs/src/components/sediment_buffer.md#push-sediment>>[3]
    frac = parcel ./ mass
    col[1,:] .+= frac .* (1.0 - bucket)
    mass -= (1.0 - bucket)
    n = floor(Int64, mass)
    
    col[n+2:end,:] .= col[1:end-n-1,:]
    # ~/~ end
    # ~/~ begin <<docs/src/components/sediment_buffer.md#push-sediment>>[4]
    na = [CartesianIndex()]
    col[2:n+1,:] .= frac[na,:]
    mass -= n
    col[1,:] .= frac .* mass
    # ~/~ end
end
# ~/~ end
# ~/~ begin <<docs/src/components/sediment_buffer.md#sediment-stack-impl>>[1]
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
    # ~/~ begin <<docs/src/components/sediment_buffer.md#pop-sediment>>[init]
    bucket = sum(col[1,:])
    @assert bucket >= 0.0
    
    if Δ < bucket
      return pop_fraction(col, Δ)
    end
    # ~/~ end
    # ~/~ begin <<docs/src/components/sediment_buffer.md#pop-sediment>>[1]
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
    # ~/~ end
end
# ~/~ end

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
# ~/~ end
