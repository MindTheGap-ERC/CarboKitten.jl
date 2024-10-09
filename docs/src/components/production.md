# Production

The `Production` module specifies the production rate following the model by Bosscher & Schlager 1992 [Bosscher1992](@cite).

``` {.julia file=src/Components/Production.jl}
@compose module Production
    @mixin WaterDepth, FaciesBase
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

``` {.julia file=src/Components/CAProduction.jl}
@compose module CAProduction
    @mixin TimeIntegration, CellularAutomaton, Production
    using ..Common
    using ..Production: production_rate
    using ..WaterDepth: water_depth

    function production(input::AbstractInput)
        w = water_depth(input)
        na = [CartesianIndex()]
        output = Array{Amount, 3}(undef, n_facies(input), input.box.grid_size...)

        w = water_depth(input)
        p(f, w) = production_rate(input.insolation, input.facies[f], w) .* input.time.Δt

        return function(state::AbstractState)
            for f = 1:n_facies(input)
                output[f, :, :] = ifelse.(state.ca .== f, p.(f, w(state)), 0.0u"m")
            end
            return output
        end 
    end
end
```
