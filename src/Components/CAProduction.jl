# ~/~ begin <<docs/src/components/production.md#src/Components/CAProduction.jl>>[init]
@compose module CAProduction
    @mixin TimeIntegration, CellularAutomaton, Production
    using ..Common
    using ..Production: capped_production, insolation
    using ..TimeIntegration: time
    using ..WaterDepth: water_depth

    function production(input::AbstractInput)
    w = water_depth(input)
    insolation_func = insolation(input)
    output_ = Array{Amount, 3}(undef, n_facies(input), input.box.grid_size...)
    n_f = n_facies(input)
    facies = input.facies
    dt = input.time.Δt

    function p(state::AbstractState, wd::AbstractMatrix)::Array{Amount,3}
        output = output_
        t = time(input.time, state)
        I = insolation_func(state)

        for i in eachindex(IndexCartesian(), wd)
            Ii = I isa AbstractArray ? I[i] : I
            active_facies = state.ca[i]
for f in 1:n_f
    output[f, i[1], i[2]] =
        f == active_facies ?
        capped_production(facies[f], Ii, wd[i], t, dt) :
        0.0u"m"
end
        end

        return output
    end

    @inline p(state::AbstractState) = p(state, w(state))
    return p
end
end

# ~/~ end