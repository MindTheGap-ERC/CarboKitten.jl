module Denudation

using Unitful

import ..BoundaryTrait: Boundary
import ..Stencil: stencil
import ..Config: Box
import ..Burgess2013.Config: Facies
import ..InputConfig: Input, DenudationType

include("./Denudation/CarbDissolution.jl")
include("./Denudation/EmpericalDenudation.jl")
include("./Denudation/PhysicalErosion.jl")

import ..Denudation.CarbDissolution: dissolution
import ..Denudation.EmpericalDenudation: emperical_denudation, slope_kernel
import ..Denudation.PhysicalErosion: physical_erosion, mass_erosion, total_mass_redistribution
export denudation, calculate_redisttribution

# configuration 

#abstract type DenudationType end

mutable struct State
    time::typeof(1.0u"Myr")
    height::Array{typeof(1.0u"m"),2}
end

@kwdef struct Dissolution <: DenudationType
    temp
    precip
    pco2
    reactionrate::Float64
end

@kwdef struct EmpericalDenudationParam <: DenudationType
    precip
end

@kwdef struct PhysicalErosionParam <: DenudationType
    erodability::Float64
end

struct NoDenudation <: DenudationType end

function water_depth(s::State)
    sea_level = input.sea_level(s.time) .* u"m"
    s.height .- sea_level
end
# generic function
"""
    denudation(param, state)

Computes the denudation for a single time-step, given denudation parameters
`param` and a simulation state `state`. `param` should have a `DenudationType`
type and `state` should contain the `height` property and `sealevel`.
"""
function denudation(::Type{BT},param::DT, water_depth::Any, slope, facies::Facies) where {BT <: Boundary, DT <: DenudationType}
end

#specific functions

function denudation(box::Box{BT},p::NoDenudation, water_depth::Any, slope, facies::Facies) where {BT <: Boundary}
    return (zeros(Float64, box.grid_size...) * u"m/kyr") 
end

function denudation(box::Box{BT}, p::EmpericalDenudationParam, water_depth::Any, slope, facies::Facies) where {BT<: Boundary} 
    return (emperical_denudation.(p.precip, slope) .* u"m/kyr") 
end

function denudation(box::Box{BT}, p::Dissolution, water_depth::Any, slope, facies::Facies) where {BT<: Boundary}
    return(dissolution(p.temp, p.precip, p.pco2, p.reactionrate, water_depth, facies) * u"m/kyr") 
end

function denudation(box::Box{BT}, p::PhysicalErosionParam, water_depth::Any, slope, facies::Facies) where {BT<: Boundary}
    # This needs transport feature to be merged so that we know the facies type of the
    # top most layer. What follows should still be regarded as pseudo-code.
    # We need to look into this further.
    denudation_amount = physical_erosion.(slope, facies.infiltration_coefficient, p.erodability)
    return (denudation_amount * u"m/kyr") 
end

#function for redistribution
function calculate_redistribution(::Type{BT},param::DT,water_depth,slope,facies) where {BT <: Boundary, DT <: DenudationType}
end

function calculate_redistribution(box::Box{BT},p::NoDenudation,water_depth,slope,facies) where {BT <: Boundary}
    return (zeros(typeof(0.0u"m/kyr"),box.grid_size...))
end

function calculate_redistribution(box::Box{BT},p::Dissolution,water_depth,slope,facies) where {BT <: Boundary}
    return (zeros(typeof(0.0u"m/kyr"),box.grid_size...))
end

function calculate_redistribution(box::Box{BT},p::EmpericalDenudationParam,water_depth,slope,facies) where {BT <: Boundary}
    @show p.precip
    return (zeros(typeof(0.0u"m/kyr"),box.grid_size...))
end

function calculate_redistribution(box::Box{BT},p::PhysicalErosionParam,water_depth, slope,facies) where {BT <: Boundary}
    redis = mass_erosion(Float64, BT, slope,(3,3),water_depth,box.phys_scale ./u"m",facies,p.erodability)
    redistribution = total_mass_redistribution(redis, slope)
    return (redistribution .* u"m/kyr")
end

"""
state needs height array
"""


end