# ~/~ begin <<docs/src/denudation/physical_erosion.md#src/Denudation/PhysicalErosionMod.jl>>[init]
module PhysicalErosionMod

import ..Abstract: DenudationType, denudation, redistribution
using ...Burgess2013: Facies
using ...Stencil: Boundary, Periodic, offset_value, offset_index, stencil
using ...BoundaryTrait
using ...Config: Box

using Unitful

@kwdef struct PhysicalErosion <: DenudationType
    erodability::typeof((1.0u"m/yr"))
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
    erodability = p.erodability ./ u"m/yr"
    denudation_amount = physical_erosion.(slope, facies.infiltration_coefficient, erodability)
    return (denudation_amount .* u"m/kyr")
end

function redistribution(box::Box{BT}, p::PhysicalErosion, water_depth, slope, inf) where {BT<:Boundary}
    erodability = p.erodability ./ u"m/yr"
    redis = mass_erosion(Float64, BT, slope, (3, 3), water_depth, box.phys_scale ./ u"m", inf, erodability)
    redistribution = total_mass_redistribution(redis, slope, BT)
    return (redistribution .* u"m/kyr")
end

end
# ~/~ end