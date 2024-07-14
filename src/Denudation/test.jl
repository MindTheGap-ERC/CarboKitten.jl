
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
		if s[i] > 0.0
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
	redis = zeros(Float64,(3,3,size(w)...))
	for i in CartesianIndices(w)
        stencil = zeros(T, n)
        for (k, Δi) in enumerate(CartesianIndices(stencil_shape))
            stencil[k] = offset_value(BT, w, i, Δi)
        end
		redis[:,:,i] .= redistribution_kernel(stencil,csz) .* physical_erosion(slope[i], inf[i], erodability)
        @show stencil redistribution_kernel(stencil,csz) physical_erosion(slope[i], inf[i], erodability)
    end
	return redis 
      					
end
#=
function add_up_kernel(redis::Array{Float64})
	r = 0.0
	r = redis[3,3,1,1] + redis[3,2,1,2] + redis[3,1,1,3]+
	redis[2,3,2,1] + redis[2,2,2,2] + redis[2,1,2,3]+
	redis[1,3,3,1] + redis[1,2,3,2] + redis[1,1,3,3]
end

function total_mass_redistribution(redis::Array{Float64},slope::Any,::Type{BT}) where {BT <: Boundary}
	s = zeros(Float64,3,3,3,3)
	shape = size(s)
	result = zeros(Float64,size(slope))
	for i in CartesianIndices(s)
		for (k, Δi) in enumerate(CartesianIndices(s))
			combined_index1 = CartesianIndex(mod1.(Tuple(i), shape)...)
            combined_index2 = CartesianIndex(mod1.(Tuple(Δi), shape)...)
			s[combined_index2,combined_index1] = offset_value(BT, redis, combined_index1, combined_index2)
		end
		result[i] = add_up_kernel(s)
	end
	return result
end
=#
function total_mass_redistribution(redis::Array{Float64},slope::Any,::Type{BT}) where {BT <: Boundary}
	mass = zeros(Float64,size(slope))
    for i in CartesianIndices(slope)
        for idx in CartesianIndices(redis)
            if offset_index(BT, size(slope), CartesianIndex(idx[3],idx[4]), CartesianIndex(idx[1]-2,idx[2]-2)) == i
			mass[i] += redis[idx]
            end
		#if idx[1] + idx[3] -1 == i[1] && idx[2] + idx[4] -1 == i[2]
		#result[i] += redis[idx]
		end
	end
	return mass
end

water_depth = -rand(5,5)
slope = rand(5,5) * 4
inf = 0.5*ones(5,5)
rs = zeros(5,5)
for i in CartesianIndices(slope)
    rs[i]=physical_erosion(slope[i],inf[i],0.23)
end
return rs
redis = mass_erosion(Float64, Periodic{2}, slope, (3,3), water_depth, 1.0,inf,0.23)
mass = total_mass_redistribution(redis,slope,Periodic{2})

