module Denudation

import ..BoundaryTrait
import ..Stencil

include("Denudation/CarbDissolution.jl")
include("Denudation/EmpericalDenudation.jl")
include("Denudation/PhysicalErosion.jl")

export denudation

# configuration 
abstract type DenudationType end

struct Dissolution <: DenudationType
    temp
    precip
    pco2
    alpha::Float64
end

struct EmpericalDenudationParam <: DenudationType
    precip
end

struct PhysicalErosionParam <: DenudationType
    erodability::Float64
end

struct NoDenudation <: DenudationType end

# generic function
"""
    denudation(param, state)

Computes the denudation for a single time-step, given denudation parameters
`param` and a simulation state `state`. `param` should have a `DenudationType`
type and `state` should contain the `height` property and `sealevel`.
"""
function denudation(box::Box{BT}, param::DT, state) where {DT <: DenudationType, BT <: BoundaryType}
end

function water_depth(state)
    state.height .- state.sea_level
end
#specific functions

function denudation(box::Box{BT}, p::EmpericalDenudationParam, state) where {BT <: BoundaryType}
    slope = zeros(Float64, box.grid_size...)
    slopefn = stencil(Float64, box, (3, 3), slope_kernel)
    w = water_depth(state)
    slopefn(w, slope, box.phys_scale) # slope is calculated with square so no need for -w
    return (emperical_denudation.(p.precip, slope), nothing)
end

function denudation(box::Box{BT}, p::Dissolution, state) where {BT <: BoundaryType}
    dissolution(p.temp, p.precip, p.co2, p.alpha, water_depth, facies)
end

function denudation(box::Box{BT}, p::PhysicalErosionParam, state) where {BT <: BoundaryType}
    # This needs transport feature to be merged so that we know the facies type of the
    # top most layer. What follows should still be regarded as pseudo-code.
    # We need to look into this further.
    slope = zeros(Float64, box.grid_size...)
    slopefn = stencil(Float64, box, (3, 3), slope_kernel)
    w = water_depth(state)
    slopefn(w, slope, box.phys_scale) # slope is calculated with square so no need for -w
    facies = (infiltration_rate = 1) # temporary place-holder
    
    redis = mass_erosion(Float64,Periodic{2},slope,(3,3),w,input.phys_scale,input.facies.inf)
    redistribution = total_mass_redistribution(redis, slope)
    denudation_amount = physical_erosion.(p.temp, p.precip, p.co2, p.alpha, w, slope, facies)

    return (denudation_amount, redistribution)
end

end