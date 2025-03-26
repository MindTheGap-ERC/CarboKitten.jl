# ~/~ begin <<docs/src/components/production.md#src/Components/CAProduction.jl>>[init]
@compose module CAProduction
    @mixin TimeIntegration, CellularAutomaton, Production
    using ..Common
    using ..Production: production_rate, insolation
    using ..WaterDepth: water_depth
    using Logging

    function production(input::AbstractInput)
        w = water_depth(input)
        na = [CartesianIndex()]
        output = Array{Amount, 3}(undef, n_facies(input), input.box.grid_size...)

        w = water_depth(input)
        s = insolation(input)
        n_f = n_facies(input)
        facies = input.facies
        Δt = input.time.Δt

        return function(state::AbstractState)
            insolation = s(state)
            for f = 1:n_f
                output[f, :, :] = ifelse.(
                    state.ca .== f,
                    production_rate.(insolation, (facies[f],), w(state)) .* Δt,
                    0.0u"m")
            end
            return output
        end
    end
end
# ~/~ end
