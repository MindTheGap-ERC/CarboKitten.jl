# ~/~ begin <<docs/src/components/denudation.md#src/Components/Denudation.jl>>[init]
@compose module Denudation
    using CarboKitten.Denudation

    @mixin Boxes, WaterDepth, SedimentBuffer, FaciesBase

    @kwdef struct Facies <: AbstractFacies
        reactive_surface::typeof(1.0u"m^2/m^3") #reactive surface
        mass_density::typeof(1.0u"kg/m^3") #density of different carb factory
        infiltration_coefficient::Float64 #infiltration coeff
        erodability::typeof(1.0u"m/yr")
    end

    @kwdef struct Input <: AbstractInput
        denudation::DenudationType
    end
end
# ~/~ end
