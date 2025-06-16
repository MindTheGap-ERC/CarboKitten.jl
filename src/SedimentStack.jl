# ~/~ begin <<docs/src/components/sediment_buffer.md#src/SedimentStack.jl>>[init]
module SedimentStack

export push_sediment!, pop_sediment!, peek_sediment

# ~/~ begin <<docs/src/components/sediment_buffer.md#sediment-stack-impl>>[init]
function push_sediment!(col::AbstractMatrix{F}, parcel::AbstractVector{F}) where F <: Real
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

    @assert mass < bucket "pop_fraction can only pop from the top cell: $(col), $(Δ)"
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

function pop_sediment!(cols::AbstractArray{F, 4}, amount::AbstractArray{F, 2}, out::AbstractArray{F, 3}) where F <: Real
  @views for i in CartesianIndices(amount)
      out[:, i[1], i[2]] = pop_sediment!(cols[:, :, i[1], i[2]], amount[i[1], i[2]])
  end
end

end # module
# ~/~ end
