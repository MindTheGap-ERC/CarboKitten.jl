# ~/~ begin <<docs/src/visualization.md#ext/ProductionCurve.jl>>[init]
module ProductionCurve

using Makie
using Unitful

import CarboKitten.Visualization: production_curve!, production_curve
using CarboKitten.Components.Production: production_rate

function production_curve!(ax, input)
    ax.title = "production at $(sprint(show, input.insolation; context=:fancy_exponent=>true))"
    ax.xlabel = "production [m/Myr]"
    ax.ylabel = "depth [m]"
    ax.yreversed = true

    for f in input.facies
        depth = (0.1:0.1:50.0)u"m"
        prod = [production_rate(input.insolation, f, d) for d in depth]
        lines!(ax, prod / u"m/Myr", depth / u"m")
    end
end

function production_curve(input)
    fig = Figure()
    ax = Axis(fig[1, 1])
    production_curve!(ax, input)
    fig
end

end
# ~/~ end
