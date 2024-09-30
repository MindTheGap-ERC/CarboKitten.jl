# ~/~ begin <<docs/src/denudation/denudation.md#src/Denudation/Abstract.jl>>[init]
module Abstract

using ...BoundaryTrait: Boundary
using ...Config: Box

using Unitful

abstract type DenudationType end

"""
    denudation(box, param, state)

FIXME Computes the denudation for a single time-step, given denudation parameters `param` and a simulation state `state`. `param` should have a `DenudationType` type and `state` should contain the `height` property and `sealevel`.
"""
function denudation(input)
    function (state, water_depth, slope)
        #denu_result = denudation(input.box, input.denudation, water_depth, slope, input.facies)
        return denudation(input.box, input.denudation, water_depth, slope, input.facies,state)
    end
end

"""
    denudation(box::Box, param::DenudationType, water_depth, slope, facies)

Computes the amount of denudation. This function is called on a pixel by pixel basis, so all arguments can be assumed to be scalar. The `param` argument should be of a subtype of `DenudationType` containing all the input parameters for this specific denudation model.
"""
function denudation(box::Box, param::DenudationType, water_depth, slope, facies)
    error("Abstract `denudation` function called.")
end

"""
    redistribution()

FIXME
"""
function redistribution(input)
        redistribution_mass::Union{Array{typeof(1.0u"m/kyr"),2}, Nothing} = nothing
    function (state, water_depth, denudation_mass)
        redi_result = redistribution(input.box, input.denudation, denudation_mass, water_depth)
        if redi_result !== nothing
            redistribution_mass = Array{typeof(1.0u"m/kyr")}(undef, size(state.ca)...)
            redistribution_mass = redi_result
        end
        return redistribution_mass
    end
end

function redistribution(box::Box, param::DenudationType, denudation_mass, water_depth)
    error("Abstract `redistribution` function called.")
end

end  # module
# ~/~ end
