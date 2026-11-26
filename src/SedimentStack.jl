# ~/~ begin <<docs/src/components/sediment_buffer.md#src/SedimentStack.jl>>[init]
module SedimentStack

export push_sediment!, pop_sediment!, peek_sediment

# ~/~ begin <<docs/src/components/sediment_buffer.md#sediment-stack-impl>>[init]
@inline function shift_down!(col::AbstractMatrix, k::Int)
    s = size(col, 1)
    @inbounds for j in axes(col, 2), i in s:-1:k+1
        col[i, j] = col[i - k, j]
    end
end

@inline function shift_up!(col::AbstractMatrix, k::Int)
    s = size(col, 1)
    @inbounds for j in axes(col, 2), i in 1:s-k
        col[i, j] = col[i + k, j]
    end
end
# ~/~ end
# ~/~ begin <<docs/src/components/sediment_buffer.md#sediment-stack-impl>>[1]
function push_sediment!(col::AbstractMatrix{F}, parcel::AbstractVector{F}) where F <: Real
    @assert size(col, 2) == length(parcel) "column $(size(col)) doesn't match parcel $(size(parcel))"
    # ~/~ begin <<docs/src/components/sediment_buffer.md#push-sediment>>[init]
    mass = sum(parcel)
    if mass == 0.0
        return
    end
    
    if mass > size(col)[1]
        @warn "pushing a very large parcel of sediment: $mass times depositional resolution"
        col .= transpose(parcel ./ mass)
        return
    end
    # ~/~ end
    # ~/~ begin <<docs/src/components/sediment_buffer.md#push-sediment>>[1]
    bucket = sum(@view col[1, :])
    @assert bucket >= 0.0 && bucket <= 1.0
    # ~/~ end
    # ~/~ begin <<docs/src/components/sediment_buffer.md#push-sediment>>[2]
    if bucket + mass < 1.0
        @views col[1, :] .+= parcel
        return
    end
    # ~/~ end
    # ~/~ begin <<docs/src/components/sediment_buffer.md#push-sediment>>[3]
    total = mass
    space = 1.0 - bucket
    @views col[1, :] .+= parcel ./ total .* space
    mass -= space
    n = floor(Int64, mass)
    
    shift_down!(col, n + 1)
    # ~/~ end
    # ~/~ begin <<docs/src/components/sediment_buffer.md#push-sediment>>[4]
    mass -= n
    for j in axes(col, 2)
        frac = parcel[j] / total
        col[2:n+1, j] .= frac
        col[1, j] = frac * mass
    end
    # ~/~ end
end
# ~/~ end
# ~/~ begin <<docs/src/components/sediment_buffer.md#sediment-stack-impl>>[2]
@inline function pop_fraction!(parcel::AbstractVector{F}, col::AbstractMatrix{F}, mass::F) where F <: Real
    bucket = sum(@view col[1, :])
    if mass == 0 || bucket == 0
        return parcel
    end

    @assert mass < bucket "pop_fraction can only pop from the top cell: $(col), $(mass)"
    r = mass / bucket
    for j in axes(col, 2)
        p = r * col[1, j]
        parcel[j] += p
        col[1, j] -= p
    end
    return parcel
end

function pop_sediment!(col::AbstractMatrix{F}, Δ::F, parcel::AbstractVector{F}) where F <: Real
    # ~/~ begin <<docs/src/components/sediment_buffer.md#pop-sediment>>[init]
    bucket = sum(@view col[1, :])
    @assert bucket >= 0.0
    
    if Δ < bucket
        parcel .= zero(F)
        return pop_fraction!(parcel, col, Δ)
    end
    # ~/~ end
    # ~/~ begin <<docs/src/components/sediment_buffer.md#pop-sediment>>[1]
    @views parcel .= col[1, :]
    Δ -= bucket
    n = floor(Int64, Δ)
    
    if n > (size(col)[1] - 2)
        @error "too much material popped of the stack: Δ = $Δ"
        for j in axes(col, 2)
            parcel[j] = sum(@view col[:, j])
        end
        col .= 0.0
        return parcel
    end
    
    for j in axes(col, 2)
        parcel[j] += sum(@view col[2:n+1, j])
    end
    shift_up!(col, n + 1)
    col[end-n:end, :] .= 0
    Δ -= n
    
    return pop_fraction!(parcel, col, Δ)
    # ~/~ end
end

pop_sediment!(col::AbstractMatrix{F}, Δ::F) where F <: Real =
    pop_sediment!(col, Δ, zeros(F, size(col, 2)))
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

  last_bucket = sum(col[n+2,:])
  if last_bucket > 0
    last_bit = (Δ / last_bucket) .* col[n+2,:]
    parcel .+= last_bit
  end

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

function pop_sediment!(cols::AbstractArray{F, 4}, amount::AbstractArray{F, 2}, out::AbstractArray{F, 3}) where F <: Real
  @views for i in CartesianIndices(amount)
      pop_sediment!(cols[:, :, i[1], i[2]], amount[i[1], i[2]], out[:, i[1], i[2]])
  end
end

end # module
# ~/~ end
