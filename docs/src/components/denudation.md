# Denudation component

``` {.julia file=src/Components/Denudation.jl}
@compose module Denudation
    using CarboKitten.Denudation

    @mixin Boxes, WaterDepth, SedimentBuffer, FaciesBase

    @kwdef struct Facies <: AbstractFacies
        reactive_surface::typeof(1.0u"m^2/m^3") #reactive surface
        mass_density::typeof(1.0u"kg/m^3") #density of different carb factory
        infiltration_coefficient::Float64 #infiltration coeff
    end

    @kwdef struct Input <: AbstractInput
        denudation::DenudationType
    end
end
```
