# Physical erosion and sediment redistribution
This method not only considers the amount of materials that have been removed, but also how the eroded materials being distributed to the neighboring regions depending on slopes on each direction.

## Physical erosion
The equations used to estimate how much material could one cell provide to the lower cells is described underneath. The equation is found in [tucker_channel-hillslope_2001](@cite). We choose this equation mainly because it specifically deals with bedrock substrates instead of loose sediments. In the equation, $k_v$ is erodibility, and the default value is 0.23 according to the paper. $(1 - I_f)$ indicates run-off generated in one cell and slope is the slope calculated based on [ArcGis: how slope works](https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/how-slope-works.htm). Note that the algorithms to calculate slope does not work on depressions.

$$D_{phys} = -k_v * (1 - I_f)^{1/3} |\nabla h|^{2/3}$$

``` {.julia #physical-erosion}
function physical_erosion(slope::Float64, Facies::facies)
    local kv = 0.23 #very arguable paramster
    #stencil(Float64,Reflected{2},(3,3),function(w)
    -kv .* (1-Facies.inf).^(1/3) .* slope.^(2/3)
end
```

## Redistribution of sediments

The redistribution of sediments after physical erosion is based on [van_de_wiel_embedding_2007](@cite): the eroded sediments calculated from the above equation are distributed to the neighboring 8 cells according to the slopes (defined as elevation differences/horizontal differences) towards each direction. The amount of sediments of one cell received is calculated by three functions below:

### Find the kernel to calculate redistibution co-efficient for the neighboring 8 cells depending on slopes

``` {.julia}
module Erosion

<<physical-erosion>>
<<erosion-transport>>
end  # module
```

``` {.julia #erosion-transport}
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
```

### Find out how much sediments would distributed to the neighboring 8 cells

``` {.julia #erosion-transport}
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
```

### How much sediment would one cell receive in total

``` {.julia #erosion-transport}
function total_mass_redistribution(redis::Array{Float64},slope::Any,::Type{BT}) where {BT <: Boundary}
	mass = zeros(Float64,size(slope))
    for i in CartesianIndices(slope)
        for idx in CartesianIndices(redis)
            if offset_index(BT, size(slope), CartesianIndex(idx[3],idx[4]), CartesianIndex(idx[1]-2,idx[2]-2)) == i
            @show i
			mass[i] += redis[idx]
            end
		end
	end
	return mass
end
```

## Current implementation

``` {.julia file=src/Denudation/PhysicalErosionMod.jl}
module PhysicalErosionMod

import ..Abstract: DenudationType, denudation, redistribution
using ...Burgess2013: Facies
using ...Stencil: Boundary, Periodic, offset_value, offset_index, stencil
using ...BoundaryTrait
using ...Config: Box

using Unitful

@kwdef struct PhysicalErosion <: DenudationType
    erodability::Float64
end

function physical_erosion(slope::Any, inf::Any, erodability::Float64)
    -1 * -erodability .* (1 - inf) .^ (1 / 3) .* slope .^ (2 / 3)
end

#erodability = 0.23

function redistribution_kernel(w::Array{Float64}, cellsize::Float64)
    s = zeros(Float64, (3, 3))
    s[1, 1] = -(w[1, 1] - w[2, 2]) / cellsize
    s[1, 2] = -(w[1, 2] - w[2, 2]) / cellsize / sqrt(2)
    s[1, 3] = -(w[1, 3] - w[2, 2]) / cellsize
    s[2, 1] = -(w[2, 1] - w[2, 2]) / cellsize / sqrt(2)
    s[2, 2] = -(w[2, 2] - w[2, 2]) / cellsize
    s[2, 3] = -(w[2, 3] - w[2, 2]) / cellsize / sqrt(2)
    s[3, 1] = -(w[3, 1] - w[2, 2]) / cellsize
    s[3, 2] = -(w[3, 2] - w[2, 2]) / cellsize / sqrt(2)
    s[3, 3] = -(w[3, 3] - w[2, 2]) / cellsize

    for i in CartesianIndices(s)
        if s[i] > 0
            continue
        else
            s[i] = 0.0
        end
    end
    sumslope = sum(s)

    if sumslope == 0.0
        zeros(Float64, (3, 3))
    else
        s ./ sumslope
    end
end

function mass_erosion(::Type{T}, ::Type{BT}, slope::Any, n::NTuple{dim,Int}, w::Array{Float64}, csz::Float64, inf::Any, erodability) where {T,dim,BT<:Boundary{dim}}
    m = n .÷ 2
    stencil_shape = range.(.-m, m)
    stencil = zeros(T, n)
    redis = zeros(Float64, (3, 3, size(w)...))
    for i in CartesianIndices(w)
        for (k, Δi) in enumerate(CartesianIndices(stencil_shape))
            stencil[k] = offset_value(BT, w, i, Δi)
        end
        redis[:, :, i] .= redistribution_kernel(stencil, csz) .* physical_erosion(slope[i], inf[i], erodability)
    end
    return redis
end

function total_mass_redistribution(redis::Array{Float64}, slope::Any, ::Type{BT}) where {BT<:Boundary}
    mass = zeros(Float64, size(slope))
    for i in CartesianIndices(slope)
        for idx in CartesianIndices(redis)
            if offset_index(BT, size(slope), CartesianIndex(idx[3], idx[4]), CartesianIndex(idx[1] - 2, idx[2] - 2)) == i
                mass[i] += redis[idx]
            end
            #if idx[1] + idx[3] -1 == i[1] && idx[2] + idx[4] -1 == i[2]
            #result[i] += redis[idx]
        end
    end
    return mass
end

function denudation(::Box, p::PhysicalErosion, water_depth::Any, slope, facies)
    # This needs transport feature to be merged so that we know the facies type of the
    # top most layer. What follows should still be regarded as pseudo-code.
    # We need to look into this further.
    denudation_amount = physical_erosion.(slope, facies.infiltration_coefficient, p.erodability)
    return (denudation_amount * u"m/kyr")
end

function redistribution(box::Box{BT}, p::PhysicalErosion, water_depth, slope, inf) where {BT<:Boundary}
    redis = mass_erosion(Float64, BT, slope, (3, 3), water_depth, box.phys_scale ./ u"m", inf, p.erodability)
    redistribution = total_mass_redistribution(redis, slope, BT)
    return (redistribution .* u"m/kyr")
end

end
```