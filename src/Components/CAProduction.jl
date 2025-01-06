# ~/~ begin <<docs/src/components/production.md#src/Components/CAProduction.jl>>[init]
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
# ~/~ end
