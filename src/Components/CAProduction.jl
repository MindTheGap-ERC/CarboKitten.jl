# ~/~ begin <<docs/src/components/production.md#src/Components/CAProduction.jl>>[init]
@compose module CAProduction
    @mixin TimeIntegration, CellularAutomaton, Production
    using ..Common
    using ..Production: insolation
    using ..WaterDepth: water_depth
    using ...Production: production_profile, capped_production
    using Logging

    function production(input::AbstractInput)
        w = water_depth(input)
        na = [CartesianIndex()]
        output_ = Array{Amount, 3}(undef, n_facies(input), input.box.grid_size...)

        w = water_depth(input)
        s = insolation(input)
        n_f = n_facies(input)
        facies = input.facies
        active_facies = [i for (i, f) in pairs(facies) if f.active]
        global_facies = [i for (i, f) in pairs(facies) if !f.active]
        dt = input.time.Δt
        # Having this a Tuple should make things type stable?
        production_specs = ((production_profile(input, f.production) for f in facies)...,)

        function p(state::AbstractState, wd::AbstractMatrix)::Array{Amount,3}
            output::Array{Amount, 3} = output_
            insolation::typeof(1.0u"W/m^2") = s(state)
            for i in eachindex(IndexCartesian(), wd)
                for f in eachindex(facies)
                    if facies[f].active
                        output[f, i[1], i[2]] = f != state.ca[i] ? 0.0u"m" :
                            capped_production(production_specs[f], insolation, wd[i], dt)
                    else
                        output[f, i[1], i[2]] =
                            capped_production(production_specs[f], insolation, wd[i], dt)
                    end
                end
            end
            return output
        end

        @inline p(state::AbstractState) = p(state, w(state))

        return p
    end
end
# ~/~ end
