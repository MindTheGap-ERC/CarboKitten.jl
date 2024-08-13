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

    s[s .< 0.0] .= 0.0
    sumslope = sum(s)

    if sumslope == 0.0
        zeros(Float64, (3, 3))
    else
        s ./ sumslope
    end
end

function mass_erosion(box::Box{BT,dim}, denudation_mass, water_depth::Array{Float64}, i::CartesianIndex) where {BT<:Boundary{2}, dim}
    wd = zeros(Float64, 3, 3)
    for (k, Δi) in enumerate(CartesianIndices((-1:1, -1:1)))
        wd[k] = offset_value(BT, water_depth, i, Δi)
    end
    return redistribution_kernel(wd, box.phys_scale |> in_units_of(u"m")) .* denudation_mass[i]
end

function total_mass_redistribution(box, denudation_mass, water_depth)
    mass = zeros(Float64, box.grid_size...)
    for i in CartesianIndices(mass)
        redis = mass_erosion(box, denudation_mass, water_depth, i)
        for subidx in CartesianIndices((-1:1, -1:1))
            target = offset_index(BT, size(slope), i, subidx)
            mass[target] += redis[2+subidx[1], 2+subidx[2]]
        end
    end
    return mass
end

function denudation(::Box, p::PhysicalErosion, water_depth::Any, slope, facies)
    # This needs transport feature to be merged so that we know the facies type of the
    # top most layer. What follows should still be regarded as pseudo-code.
    # We need to look into this further.
    erodability = p.erodability ./ u"m/yr"
    denudation_mass = physical_erosion.(slope, facies.infiltration_coefficient, erodability)
    return (denudation_mass .* u"m/kyr")
end

function redistribution(input)
    function (state)
    end
end

function redistribution(box::Box{BT}, ::PhysicalErosion, denudation_mass, water_depth) where {BT<:Boundary}
    return total_mass_redistribution(box, denudation_mass, water_depth) * u"m/kyr"
end

end
# ~/~ end
