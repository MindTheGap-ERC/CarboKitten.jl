# ~/~ begin <<docs/src/Denudation.md#src/Erosion.jl>>[init]
module Erosion

# ~/~ begin <<docs/src/Denudation.md#physical-erosion>>[init]
function physical_erosion(slope::Float64, Facies::facies)
    local kv = 0.23 #very arguable paramster
    #stencil(Float64,Reflected{2},(3,3),function(w)
    -kv .* (1-Facies.inf).^(1/3) .* slope.^(2/3)
end
# ~/~ end
# ~/~ begin <<docs/src/Denudation.md#erosion-transport>>[init]
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
# ~/~ end
# ~/~ begin <<docs/src/Denudation.md#erosion-transport>>[1]
function mass_erosion(::Type{T},::Type{BT},slope::Matrix{Float64},n::NTuple{dim,Int}) where {T, dim, BT <: Boundary{dim}}
	m = n .÷ 2
    stencil_shape = range.(.-m, m)
    stencil = zeros(T, n)
	redis = zeros(Float64,(3,3,size(slope)...))
	local inf = 0.5
	for i in CartesianIndices(slope)
	     #println(i)
        for (k, Δi) in enumerate(CartesianIndices(stencil_shape))
			#println(Δi)
            stencil[k] = offset_value(BT, w, i, Δi)
			#println(k)
			redis[:,:,i] .= -1 .* redistribution_kernel(stencil,csz) .* physical_erosion(slope[i],inf)
        end
    end
	return redis		

end
# ~/~ end
# ~/~ begin <<docs/src/Denudation.md#erosion-transport>>[2]
function total_mass_redistribution(redis::Array{Float64},slope::Matrix{Float64})
	result = zeros(Float64,size(slope))
	for idx in CartesianIndices(redis)
		for i in CartesianIndices(slope)
			#println(idx)
			#println(i)
			#println(idx[1],idx[2],idx[3],idx[4],i[1],i[2])
		if idx[1] + idx[3] -1 == i[1] && idx[2] + idx[4] -1 == i[2]
		result[i] += redis[idx]
		end
		end
	end
	return result
end
# ~/~ end
end  # module
# ~/~ end