# ~/~ begin <<docs/src/denudation/denudation.md#src/Denudation/NoDenudationMod.jl>>[init]
"""
    module NoDenudation

Doesn't do any denudation: used for testing purposes.
"""
module NoDenudationMod

import ..Abstract: DenudationType, denudation, redistribution
using ...Boxes: Box
using Unitful

struct NoDenudation <: DenudationType end

function denudation(box::Box, p::NoDenudation, water_depth::Any, slope, facies, state)
    return nothing
end

function redistribution(box::Box, p::NoDenudation, denudation_mass, water_depth)
    return nothing
end

end
# ~/~ end