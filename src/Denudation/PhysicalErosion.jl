module PhysicalErosion

using CarboKitten.Burgess2013: Facies
using CarboKitten.Stencil: Boundary, Periodic, offset_value, offset_index, stencil
using CarboKitten.BoundaryTrait
import CarboKitten.Config: Box
export physical_erosion, mass_erosion, total_mass_redistribution

function physical_erosion(slope::Any, inf::Any, erodability::Float64)
    -1 * -erodability .* (1-inf).^(1/3) .* slope.^(2/3)
end

#erodability = 0.23

function redistribution_kernel(w::Array{Float64},cellsize::Float64)
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

function mass_erosion(::Type{T},::Type{BT},slope::Any,n::NTuple{dim,Int},w::Array{Float64},csz::Float64,inf::Any,erodability) where {T, dim, BT <: Boundary{dim}}
	m = n .÷ 2
    stencil_shape = range.(.-m, m)
    stencil = zeros(T, n)
	redis = zeros(Float64,(3,3,size(w)...))
	for i in CartesianIndices(w)
        for (k, Δi) in enumerate(CartesianIndices(stencil_shape))
            stencil[k] = offset_value(BT, w, i, Δi)
			redis[:,:,i] .= redistribution_kernel(stencil,csz) .* physical_erosion(slope[i], inf[i], erodability)
        end
    end
	return redis		
				
end
	
function total_mass_redistribution(redis::Array{Float64},slope::Any)
	result = zeros(Float64,size(slope))
	for idx in CartesianIndices(redis)
		for i in CartesianIndices(slope)
		if idx[1] + idx[3] -1 == i[1] && idx[2] + idx[4] -1 == i[2]
		result[i] += redis[idx]
		end
		end
	end
	return result
end

end