"""
    module NoDenudation

Doesn't do any denudation: used for testing purposes.
"""
module NoDenudationMod

import ..Abstract: DenudationType, denudation, redistribution
using Unitful
using ...Config: Box

struct NoDenudation <: DenudationType end

function denudation(box::Box, p::NoDenudation, water_depth::Any, slope, facies)
    return (zeros(typeof(0.0u"m/kyr"), box.grid_size...))
end

function redistribution(box::Box, p::NoDenudation, water_depth, slope, facies)
    return (zeros(typeof(0.0u"m/kyr"), box.grid_size...))
end

end

