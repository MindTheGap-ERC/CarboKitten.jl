### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ d0bb3cb2-df98-488a-9c99-1d530d846e01
using Pkg; Pkg.activate("..")

# ╔═╡ 64d27eb9-e43a-4552-bdf4-3ceaac66b2aa
using CarboKitten.EmpericalDenudation

# ╔═╡ 5badffc6-a814-422c-9872-751de617fdf8
using CarboKitten

# ╔═╡ e3c85e34-0a09-4ae8-a8b2-abc93862c835
using CarboKitten.Stencil

# ╔═╡ f65a48f0-cb23-11ee-379b-c13d507197da
w = rand(10,10);

# ╔═╡ 1d4bb849-b1af-4eb0-8a9f-54a01d3c7d55
function redistribution_kernel(w::Matrix{Float64},cellsize::Float64)
    s = zeros(Float64,(3,3))
	s[1,1] = -(w[1,1] - w[2,2]) / cellsize
    s[1,2] = -(w[1,2] - w[2,2]) / cellsize / sqrt(2)
    s[1,3] = -(w[1,3] - w[2,2]) / cellsize
    s[2,1] = -(w[2,1] - w[2,2]) / cellsize / sqrt(2) 
    s[2,2] = -(w[2,2] - w[2,2]) / cellsize
    s[2,3] = -(w[2,3] - w[2,2]) / cellsize / sqrt(2)
    s[3,1] = -(w[3,1] - w[2,2]) / cellsize 
    s[3,2] = -(w[3,2] - w[2,2]) / cellsize / sqrt(2)
    s[3,3] = -(w[3,3] - w[2,2]) / cellsize
	
	for i in CartesianIndices(s)
		if s[i] > 0
		   continue
		else
		   s[i] = 0.0
		end
	end
	sumslope = sum(s)
	
	if sumslope == 0.0
	zeros(Float64,(3,3))
	else
	s./sumslope
	end
end


# ╔═╡ de2f1418-984c-498f-bcf0-ecc647dbe39e
csz = 1.0

# ╔═╡ 54d5bc80-8f61-4d64-961b-cc51fe86d207
redistribution_kernel(w,csz)

# ╔═╡ d1115b0e-6f9f-4fe2-9a5d-c73c2294bb68
const kv = 0.23;

# ╔═╡ bafbe0cb-dd89-4dcf-8677-282e093998f5
function physical_erosion(slope::Float64, inf::Float64)
    -kv .* (1-inf).^(1/3) .* slope.^(2/3)
end

# ╔═╡ 597519d2-3fef-4d2c-9575-9344b26b5bbb


# ╔═╡ 07e5753e-cf2d-4d6e-aacc-cb218c4a4f43
slopefn = stencil(Float64, Periodic{2}, (3, 3), slope_kernel) #function 

# ╔═╡ e5cc2e5a-b9b7-4d04-b9f2-071ba57dce9e
slope = zeros(size(w)...);#precallocate memory

# ╔═╡ 7c5c6103-3426-42d2-94aa-56aa7ee18263
slopefn(w,slope,csz); #instance

# ╔═╡ 90cb0c58-070c-4c78-8c4a-d73506fe0374
D = ones(size(w)...);

# ╔═╡ da9201ec-6b0d-48cb-97d9-7ca6f60d169b
for i in CartesianIndices(slope)
	D[i] = physical_erosion(slope[i],0.5)
end

# ╔═╡ 80e44d1d-f2c2-414a-b3fc-d3c59fbaf94c
D[5,3]

# ╔═╡ b4dd2f71-2cf9-4483-be31-837de6de93b0
function offset_index(::Type{Periodic{dim}}, shape::NTuple{dim,Int}, i::CartesianIndex, Δi::CartesianIndex) where {dim}
    CartesianIndex(mod1.(Tuple(i + Δi), shape)...)
end

# ╔═╡ 3194a91e-9f22-4eb5-a5e6-c1b4c89484a5
function offset_value(BT::Type{B}, z::AbstractArray, i::CartesianIndex, Δi::CartesianIndex) where {dim, B <: Boundary{dim}}
    z[offset_index(BT, size(z), i, Δi)]
end

# ╔═╡ 953b7f2f-1d81-407c-919b-dc8f1eae75bd
function mass_erosion(::Type{T},::Type{BT},slope::Matrix{Float64},n::NTuple{dim,Int}) where {T, dim, BT <: Boundary{dim}}
	m = n .÷ 2
    stencil_shape = range.(.-m, m)
    stencil = zeros(T, n)
	redis = zeros(Float64,(3,3,size(slope)...))
	result = zeros(Float64,size(slope))
	local inf = 0.5
	for i in CartesianIndices(slope)
	     #println(i)
        for (k, Δi) in enumerate(CartesianIndices(stencil_shape))
			#println(Δi)
            stencil[k] = offset_value(BT, w, i, Δi)
			#println(k)
			redis[:,:,i] .= redistribution_kernel(stencil,csz) .* physical_erosion(slope[i],0.5)
        end
    end
	return redis		
				
end
	

# ╔═╡ 0a8238d9-4413-40de-81fe-9081826479b4


# ╔═╡ a1a4ad8e-1048-49c4-b627-7f205602f97c
redis = mass_erosion(Float64,Periodic{2},slope,(3,3))

# ╔═╡ e40726d5-08e1-4011-9003-46f4e80234d0
function total_mass(redis::Array{Float64},slope::Matrix{Float64})
	result = zeros(Float64,size(slope))
	for idx in CartesianIndices(redis)
		for i in CartesianIndices(slope)
			println(idx)
			println(i)
			println(idx[1],idx[2],idx[3],idx[4],i[1],i[2])
		if idx[1] + idx[3] -1 == i[1] && idx[2] + idx[4] -1 == i[2]
		result[i] += redis[idx]
		end
		end
	end
	return result
end

# ╔═╡ d55d4351-6a58-4df9-9624-d5546bfdb6bd
total_mass(redis,slope)

# ╔═╡ 7ba510f9-c5c1-4195-b7df-4f22d0e74245
function mass_add(w::Matrix{Float64},redis::Array{Float64})
	result = zeros(Float64,size(w))
	for i in 2:size(w)
		println(i)
			result[i] = redis[3,3,(Tuple(i) .+ (-1,-1))] .+ redis[3,2,(Tuple(i) .+ (-1,0))] .+ redis[3,1,(Tuple(i) .+ (-1,1))] .+ redis[2,3,(Tuple(i) .+ (0,-1))] .+ redis[2,1,(Tuple(i) .+ (0,1))] .+ redis[1,3,(Tuple(i) .+ (1,-1))] .+ redis[2,3,(Tuple(i) .+ (0,-1))] .+ redis[1,1,(Tuple(i) .+ (1,1))]
	end
end

# ╔═╡ Cell order:
# ╠═f65a48f0-cb23-11ee-379b-c13d507197da
# ╠═1d4bb849-b1af-4eb0-8a9f-54a01d3c7d55
# ╠═de2f1418-984c-498f-bcf0-ecc647dbe39e
# ╠═54d5bc80-8f61-4d64-961b-cc51fe86d207
# ╠═d1115b0e-6f9f-4fe2-9a5d-c73c2294bb68
# ╠═bafbe0cb-dd89-4dcf-8677-282e093998f5
# ╠═597519d2-3fef-4d2c-9575-9344b26b5bbb
# ╠═64d27eb9-e43a-4552-bdf4-3ceaac66b2aa
# ╠═07e5753e-cf2d-4d6e-aacc-cb218c4a4f43
# ╠═e5cc2e5a-b9b7-4d04-b9f2-071ba57dce9e
# ╠═7c5c6103-3426-42d2-94aa-56aa7ee18263
# ╠═90cb0c58-070c-4c78-8c4a-d73506fe0374
# ╠═da9201ec-6b0d-48cb-97d9-7ca6f60d169b
# ╠═80e44d1d-f2c2-414a-b3fc-d3c59fbaf94c
# ╠═d0bb3cb2-df98-488a-9c99-1d530d846e01
# ╠═5badffc6-a814-422c-9872-751de617fdf8
# ╠═e3c85e34-0a09-4ae8-a8b2-abc93862c835
# ╠═b4dd2f71-2cf9-4483-be31-837de6de93b0
# ╠═3194a91e-9f22-4eb5-a5e6-c1b4c89484a5
# ╠═953b7f2f-1d81-407c-919b-dc8f1eae75bd
# ╠═0a8238d9-4413-40de-81fe-9081826479b4
# ╠═a1a4ad8e-1048-49c4-b627-7f205602f97c
# ╠═e40726d5-08e1-4011-9003-46f4e80234d0
# ╠═d55d4351-6a58-4df9-9624-d5546bfdb6bd
# ╠═7ba510f9-c5c1-4195-b7df-4f22d0e74245
