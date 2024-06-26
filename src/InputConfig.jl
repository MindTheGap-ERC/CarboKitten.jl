module InputConfig

using Unitful
import ..Config: Box, TimeProperties
import ..Burgess2013.Config: Facies


export Input, DenudationType

abstract type DenudationType end

@kwdef struct Input
    box :: Box
    time :: TimeProperties

    sea_level       # Myr -> m
    subsidence_rate::typeof(1.0u"m/Myr")
    initial_depth  # m -> m

    facies::Vector{Facies}
    insolation::typeof(1.0u"W/m^2")

    denudationparam::DenudationType
end

end