
``` {.julia file=src/Components/Production.jl}
@compose module UniformProduction
    @mixin WaterDepth, FaciesConcept
    using ..Common
    using ..WaterDepth: water_depth

    export production_rate, uniform_production

    @kwdef struct Facies <: AbstractFacies
        maximum_growth_rate::Rate
        extinction_coefficient::typeof(1.0u"m^-1")
        saturation_intensity::Intensity
    end

    @kwdef struct Input <: AbstractInput
        insolation::Intensity
    end

    function production_rate(insolation, facies, water_depth)
        gₘ = facies.maximum_growth_rate
        I = insolation / facies.saturation_intensity
        x = water_depth * facies.extinction_coefficient
        return water_depth > 0.0u"m" ? gₘ * tanh(I * exp(-x)) : 0.0u"m/Myr"
    end

    function uniform_production(input::AbstractInput)
        w = water_depth(input)
        na = [CartesianIndex()]

        return function(state::AbstractState)
        return production_rate.(
            input.insolation,
            input.facies[:,na,na],
            w(state)[na,:,:])
        end
    end
end
```
