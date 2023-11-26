# ~/~ begin <<docs/src/submarine-transport.md#src/SedimentStack.jl>>[init]
module SedimentStack

export push_sediment!, pop_sediment!

function push_sediment!(col::AbstractMatrix{F}, parcel::AbstractVector{F}) where F <: Real
  Δ = sum(parcel)
  bucket = sum(col[1, :])

  if bucket + Δ < 1.0
    col[1,:] .+= parcel
    return
  end

  frac = parcel ./ Δ
  col[1,:] .+= frac .* (1.0 - bucket)
  Δ -= (1.0 - bucket)
  n = floor(Int64, Δ)

  col[n+2:end,:] = col[1:end-n-1,:]
  na = [CartesianIndex()]
  col[2:n+1,:] .= frac[na,:]
  Δ -= n

  col[1,:] .= frac .* Δ
end

@inline function pop_fraction(col::AbstractMatrix{F}, Δ::F) where F <: Real
  bucket = sum(col[1,:])
  @assert Δ < bucket "pop_fraction can only pop from the top cell"
  parcel = (Δ / bucket) .* col[1,:]
  col[1,:] .-= parcel
  return parcel
end

function pop_sediment!(col::AbstractMatrix{F}, Δ::F) where F <: Real  # -> Vector{F}
  bucket = sum(col[1,:])
  if Δ < bucket
    return pop_fraction(col, Δ)
  end

  parcel = copy(col[1,:])
  Δ -= bucket
  n = floor(Int64, Δ)

  parcel .+= sum(col[2:n+1,:]; dims=1)'
  col[1:end-n-1, :] = col[n+2:end, :]
  col[end-n-1:end, :] .= 0
  Δ -= n

  parcel .+= pop_fraction(col, Δ)
  return parcel
end

end # module
# ~/~ end