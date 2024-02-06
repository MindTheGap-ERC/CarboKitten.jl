# ~/~ begin <<docs/src/ca-with-production.md#src/ForwardModel.jl>>[init]
module ForwardModel

function run(state, propagator, updater)
    Channel{Frame}() do ch
        while true
            Δ = p(state)
            put!(ch, Δ)
            u(state, Δ)
        end
    end
end

end
# ~/~ end