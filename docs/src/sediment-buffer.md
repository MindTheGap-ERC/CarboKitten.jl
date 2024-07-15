# Sediment Buffers

In some models of transport it is important to remember the type of sediment facies some time in the past. Suppose sedimented material is loosened again due to submarine erosion processes. Particle sizes may depend on the sediment facies.

To compute such models efficiently we need to store the sedimentation history in a buffer. One approach is to remember several iterations into the past, but then there will be areas where sedimentation is slow and fast, and soon you'll find you have to store the entire simulation for retrieval.

The proposed data structure is one that stores the history discretized in depth. We loose a bit of precision on the thinnest sedimentation layers, but we can freely mix older and younger sediments when needed.

There are several choices on how to structure the sediment buffer. We can grow sediment from the bottom of the buffer, but that requires keeping an array of pointers to the bottom depth where new material is deposited.

Another choice is to keep the sea floor in the same layer of the buffer, and copy down sediment when a layer is full. We may need to implement both to see which is more efficient.

We define two functions `push_sediment!` and `pop_sediment!`. Given a $s \times n$ matrix, where $n$ is the number of facies types and $s$ is the depth of the stack, we can grow and shrink sediment. These functions are unit-free, setting $\Delta z$ to be equal to 1.

``` {.julia file=test/SedimentStackSpec.jl}
@testset "SedimentStack" begin
  using CarboKitten.SedimentStack: push_sediment!, pop_sediment!
  stack = zeros(Float64, 10, 3)
  push_sediment!(stack, [5.0, 0, 0])
  @test pop_sediment!(stack, 1.5) == [1.5, 0.0, 0.0]
  push_sediment!(stack, [0.0, 2.0, 0.0])   # (0 0.5) (0 1) (0.5 0.5) (1 0) ...
  @test pop_sediment!(stack, 2.0) == [0.25, 1.75, 0.0]
  @test pop_sediment!(stack, 1.5) == [1.25, 0.25, 0.0]
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
```

```@raw html
<details><summary>SedimentStack impl</summary>
```

``` {.julia file=src/SedimentStack.jl}
module SedimentStack

export push_sediment!, pop_sediment!, peek_sediment

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

function push_sediment!(sediment::AbstractArray{F, 4}, p::AbstractArray{F, 3}) where F <: Real
  _, x, y = size(p)
  @views for i in CartesianIndices((x, y))
    push_sediment!(sediment[:, :, i[1], i[2]], p[:, i[1], i[2]])
  end
end

@inline function pop_fraction(col::AbstractMatrix{F}, Δ::F) where F <: Real
  bucket = sum(col[1,:])
  @assert Δ < bucket "pop_fraction can only pop from the top cell"
  parcel = (Δ / bucket) .* col[1,:]
  col[1,:] .-= parcel
  return parcel
end

function peek_sediment(col::AbstractMatrix{F}, Δ::F) where F <: Real  # -> Vector{F}
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
```

```@raw html
</details>
```
