# ~/~ begin <<docs/src/visualization.md#ext/ProductionCurve.jl>>[init]
module ProductionCurve

import CarboKitten.Visualization: production_curve!, production_curve
using CarboKitten.Burgess2013: production_rate

function production_curve!(ax, input; loc=nothing)
end

function production_curve(input; loc=nothing)
    fig, loc = isnothing(loc) ? let fig = Figure()
        (fig, fig[1, 1])
    end : (nothing, loc)
    ax = Axis(loc, title="production at $(sprint(show, input.insolation; context=:fancy_exponent=>true))",
              xlabel="production (m/Myr)", ylabel="depth (m)", yreversed=true)
    for f in input.facies
        depth = (0.1:0.1:50.0)u"m"
        prod = [production_rate(input.insolation, f, d) for d in depth]
        lines!(ax, prod / u"m/Myr", depth / u"m")

    end
    fig
end

end
# ~/~ end