# ~/~ begin <<docs/src/components/cellular-automata.md#src/Components/CAFeedback.jl>>[init]
@compose module CAFeedback
using ..Common

@kwdef struct Facies
    minimum_production::Union{typeof(0.0u"m/Myr"),Nothing} = nothing
end

function ca_feedback(input::AbstractInput)
    production_limit = [f.minimum_production for f in input.facies]
    dt = input.time.Δt

    function (ca, production)
        for i in eachindex(IndexCartesian(), ca)
            f = ca[i]
            if f != 0 && production_limit[f] !== nothing &&
                production[f, i[1], i[2]] / dt < production_limit[f]
                ca[i] = 0
            end
        end
    end
end

end
# ~/~ end
