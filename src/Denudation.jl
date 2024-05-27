module Denudation

import ..BoundaryTrait
import ..Stencil
import ..Config: Box
import ..Burgess2013.Config: Facies
import ..InputConfig: Input, DenudationType


include("./Denudation/CarbDissolution.jl")
include("./Denudation/EmpericalDenudation.jl")
include("./Denudation/PhysicalErosion.jl")

export denudation

# configuration 

#abstract type DenudationType end
abstract type State end
struct Dissolution <: DenudationType
    temp
    precip
    pco2
    reactionrate::Float64
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
function denudation(input::Input, param::DT, state) where {DT <: DenudationType}
end

function water_depth(state::State)
    state.height .- state.sea_level
end
#specific functions

function denudation(input::Input, p::NoDenudation, state) 
    return (nothing,nothing)
end

function denudation(input::Input, p::EmpericalDenudationParam, state::State) 
    slope = zeros(Float64, input.box.grid_size...)
    slopefn = stencil(Float64, input.box{BT}, (3, 3), slope_kernel)
    w = water_depth(state)
    slopefn(w, slope, input.box.phys_scale) # slope is calculated with square so no need for -w
    return (emperical_denudation.(p.precip, slope), nothing)
end

function denudation(input::Input, p::Dissolution, state::State) 
    w = water_depth(state)
    return(dissolution(p.temp, p.precip, p.co2, p.reactionrate, w, input.facies.infiltration_coefficient),nothing)
end

function denudation(input::Input, p::PhysicalErosionParam, state::State) 
    # This needs transport feature to be merged so that we know the facies type of the
    # top most layer. What follows should still be regarded as pseudo-code.
    # We need to look into this further.
    slope = zeros(Float64, input.box.grid_size...)
    slopefn = stencil(Float64, input.box{BT}, (3, 3), slope_kernel)
    w = water_depth(state)
    slopefn(w, slope, input.box.phys_scale) # slope is calculated with square so no need for -w
    denudation_amount = physical_erosion.(slope, input.facies.infiltration_coefficient, p.erodability)
    redis = mass_erosion(Float64,input.box{BT}, slope,(3,3),w,box.phys_scale,input.facies.infiltration_coefficient,input.erodability)
    redistribution = total_mass_redistribution(redis, slope)

    return (denudation_amount, redistribution)
end


"""
state needs height array
"""


end