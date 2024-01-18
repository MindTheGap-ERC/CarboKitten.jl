### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ b9c5c039-ca94-431e-946b-a027bfc437ff
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

# ╔═╡ d3966730-5f3c-4baf-8a25-74e85436d57f
function push_sediment!(sediment::AbstractArray{F, 4}, p::AbstractArray{F, 3}) where F <: Real
  _, x, y = size(p)
  @views for i in CartesianIndices((x, y))
    push_sediment!(sediment[:, :, i[1], i[2]], p[:, i[1], i[2]])
  end
end

# ╔═╡ 56b665d6-bcee-40de-bd83-458521394144
@inline function pop_fraction(col::AbstractMatrix{F}, Δ::F) where F <: Real
  bucket = sum(col[1,:])
  @assert Δ < bucket "pop_fraction can only pop from the top cell"
  parcel = (Δ / bucket) .* col[1,:]
  col[1,:] .-= parcel
  return parcel
end

# ╔═╡ f12ddf24-b2e4-4d29-93ae-4df5e4b97e07
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

# ╔═╡ 57c7b95b-c9ea-4748-b4c8-84522101bf5c
function peek_sediment(sediment::AbstractArray{F,4}, Δ::F) where F <: Real
  _, f, x, y = size(sediment)
  out = Array{F, 3}(undef, f, x, y)
  for i in CartesianIndices((x, y))
    out[:, i[1], i[2]] = peek_sediment(@view(sediment[:, :, i[1], i[2]]), Δ)
  end
  return out
end

# ╔═╡ df0f6732-7c7b-42b0-8ae4-7693360b2177
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

# ╔═╡ 691028eb-af16-4c34-bc47-f4a746aa9461
x = let
	x = zeros(10, 3, 5, 5)
	for i in 1:10
		p = rand(3, 5, 5)
		push_sediment!(x, p)
	end
	x
end

# ╔═╡ a08024c0-cbb6-4fd6-9704-4e5c05ea79dc
sum(peek_sediment(x, 1.0); dims=1)

# ╔═╡ 3d4bc851-19be-4ae4-9039-d6d75ea992e6


# ╔═╡ a4a364fd-8ca7-450f-9589-0635edcad3cd


# ╔═╡ Cell order:
# ╠═b9c5c039-ca94-431e-946b-a027bfc437ff
# ╠═d3966730-5f3c-4baf-8a25-74e85436d57f
# ╠═56b665d6-bcee-40de-bd83-458521394144
# ╠═f12ddf24-b2e4-4d29-93ae-4df5e4b97e07
# ╠═57c7b95b-c9ea-4748-b4c8-84522101bf5c
# ╠═df0f6732-7c7b-42b0-8ae4-7693360b2177
# ╠═691028eb-af16-4c34-bc47-f4a746aa9461
# ╠═a08024c0-cbb6-4fd6-9704-4e5c05ea79dc
# ╠═3d4bc851-19be-4ae4-9039-d6d75ea992e6
# ╠═a4a364fd-8ca7-450f-9589-0635edcad3cd
