module Denudation

include("Denudation/CarbDissolution.jl")
include("Denudation/EmpericalDenudation.jl")
include("Denudation/PhysicalErosion.jl")

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
function denudation(::DT, water_depth::Float64, facies::Facies) where {DT <: DenudationType}
end

#specific functions
function denudation(p::EmpericalDenudationParam, water_depth::Float64, facies::Facies)
    s
    emperical_denudation(p.precip, water_depth, facies)
end

function denudation(p::Dissolution, water_depth::Float64, facies::Facies)
    dissolution(p.temp, p.precip, p.co2, p.alpha, water_depth, facies)
end

function denudation(p::PhysicalErosionParam, water_depth::Float64, slope::Float64, facies::Facies)
    physicalerosion(p.temp, p.precip, p.co2, p.alpha, water_depth, slope, facies)
end

end