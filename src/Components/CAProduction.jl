# ~/~ begin <<docs/src/components/production.md#src/Components/CAProduction.jl>>[init]
@compose module CAProduction
    @mixin TimeIntegration, CellularAutomaton, Production
    using ..Common
    using ..Production: production_rate
    using ..WaterDepth: water_depth

    function production(input::AbstractInput)
        w = water_depth(input)
        na = [CartesianIndex()]
        output = zeros(Amount, n_facies(input), input.box.grid_size...)

        return function(state::AbstractState)
            output[state.ca,:,:] = production_rate.(
                input.insolation,
                input.facies[state.ca],
                w(state))
            return output
        end 
    end
end
# ~/~ end
