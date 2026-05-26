# Copying Stack

We define two functions `push_sediment!` and `pop_sediment!`. Given a $s \times n$ matrix, where $n$ is the number of facies types and $s$ is the depth of the stack, we can grow and shrink sediment. These functions are unit-free, setting $\Delta z$ to be equal to 1.

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

Implementation
--------------

``` {.julia file=src/Algorithms/CopyBuffers.jl}
module CopyBuffers

module Internal

<<sediment-stack-impl>>

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

end  # module Internal

import ...Interfaces.SedimentBuffers: push_sediment!, pop_sediment!, peek_sediment
using ...Interfaces.Chunks

struct CopyBuffer{F}
    data::Array{Float64, 4}
end

function push_sediment!(buffer::CopyBuffer{F}, sediment::AbstractArray{F, 3}, chunk::Serial) where F
    @views for i in CartesianIndices(normalize(chunk.slice, buffer.data)...)
        Internal.push_sediment!(buffer.data[:, :, i[1], i[2]], p[:, i[1], i[2]])
    end
end

function pop_sediment!(buffer::CopyBuffer, amount::AbstractArray{F, 2}, out::AbstractArray{F, 3}, chunk::Serial) where F
    @views for i in CartesianIndices(normalize(chunk.slice, buffer.data)...)
        out[:, i[1], i[2]] = Internal.pop_sediment!(buffer.data[:, :, i[1], i[2]], amount[i])
    end
end

function peek_sediment(buffer::CopyBuffer, amount::AbstractArray{F, 2}, out::AbstractArray{F, 3}, chunk::Serial) where F
    @views for i in CartesianIndices(normalize(chunk.slice, buffer.data)...)
        out[:, i[1], i[2]] = Internal.peek_sediment(buffer.data[:, :, i[1], i[2]], amount[i])
    end
end

end  # module CopyBuffers
```
