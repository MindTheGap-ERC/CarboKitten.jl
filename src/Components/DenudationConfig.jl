# ~/~ begin <<docs/src/components/denudation.md#src/Components/Denudation.jl>>[init]
@compose module DenudationConfig
    @mixin Boxes, WaterDepth, SedimentBuffer, FaciesBase
    using ..Common
    using CarboKitten.Denudation: DenudationType
    using CarboKitten.BoundaryTrait

    export slope_function
    
    @kwdef struct Facies <: AbstractFacies
        reactive_surface::typeof(1.0u"m^2/m^3") #reactive surface
        mass_density::typeof(1.0u"kg/m^3") #density of different carb factory
        infiltration_coefficient::Float64 #infiltration coeff
        erodability::typeof(1.0u"m/yr")
    end

    @kwdef struct Input <: AbstractInput
        denudation::DenudationType
    end

    function slope_function(input::Input,box::Box{BT}) where {BT<:Boundary}
    return slopefn = stencil(Float64, BT, (3, 3), slope_kernel) 
    end
end
# ~/~ end
