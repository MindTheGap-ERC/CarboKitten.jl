# ~/~ begin <<docs/src/denudation/denudation.md#src/Denudation/Abstract.jl>>[init]
module Abstract

using ...BoundaryTrait: Boundary
using ...Boxes: Box

using Unitful

abstract type DenudationType end

"""
    denudation(box, param, state)

FIXME Computes the denudation for a single time-step, given denudation parameters `param` and a simulation state `state`. `param` should have a `DenudationType` type and `state` should contain the `height` property and `sealevel`.
"""
function denudation(input)
    function (state, water_depth, slope)
        return denudation(input.box, input.denudation, water_depth, slope, input.facies,state)
    end
end

"""
    denudation(box::Box, param::DenudationType, water_depth, slope, facies)

Computes the amount of denudation. This function is called on a pixel by pixel basis, so all arguments can be assumed to be scalar. The `param` argument should be of a subtype of `DenudationType` containing all the input parameters for this specific denudation model.
"""
function denudation(box::Box, param::DenudationType, water_depth, slope, facies, state)
    error("Abstract `denudation` function called.")
end

"""
    redistribution()

FIXME
"""
function redistribution(input)
    function (state, water_depth, denudation_mass)
        return redistribution(input.box, input.denudation, denudation_mass, water_depth)
    end
end

function redistribution(box::Box, param::DenudationType, denudation_mass, water_depth)
    error("Abstract `redistribution` function called.")
end

end  # module
# ~/~ end