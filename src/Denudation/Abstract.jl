# ~/~ begin <<docs/src/denudation/denudation.md#src/Denudation/Abstract.jl>>[init]
module Abstract

using ...BoundaryTrait: Boundary

abstract type DenudationType end

"""
    denudation(box, param, state)

FIXME Computes the denudation for a single time-step, given denudation parameters `param` and a simulation state `state`. `param` should have a `DenudationType` type and `state` should contain the `height` property and `sealevel`.
"""
function denudation(input)
    denudation_mass::Array{typeof(1.0u"m/kyr"),2} = Array(undef, input.box.grid_size...)

    function (state, water_depth, slope)
        for idx in CartesianIndices(state.ca)
            f = state.ca[idx]
            if f == 0
                continue
            end

            if water_depth[idx] >= 0
                (denudation_mass[idx]) = denudation(input.box, input.denudation, water_depth[idx], slope[idx], input.facies[f])
            end
        end
        return denudation_mass
    end
end

function denudation(box::Box, param::DenudationType, water_depth, slope, facies)
    error("Abstract `denudation` function called.")
end

"""
    redistribution()

FIXME
"""
function redistribution(input)
    n_facies = length(input.facies)
    redistribution_mass::Array{typeof(1.0u"m/kyr"),3} = Array(undef, n_facies, input.box.grid_size...)

    function (state, denudation_mass, slope, inf_map)
        return redistribution_mass
    end
end

function redistribution(box::Box, param::DenudationType, water_depth, slope, inf_map)
    error("Abstract `redistribution` function called.")
end

end  # module
# ~/~ end
