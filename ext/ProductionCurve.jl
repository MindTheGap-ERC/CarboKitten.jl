# ~/~ begin <<docs/src/visualization.md#ext/ProductionCurve.jl>>[init]
module ProductionCurve

using Makie
using Unitful
using HDF5

import CarboKitten.Visualization: production_curve!, production_curve
using CarboKitten.Components.Production: Facies, production_rate

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

function production_curve!(ax, g::HDF5.Group)
    ax.title = "production at $(sprint(show, input.insolation; context=:fancy_exponent=>true))"
    ax.xlabel = "production [m/Myr]"
    ax.ylabel = "depth [m]"
    ax.yreversed = true

    a = attributes(g)
    insolation = a["insolation"][]

    for i in 1:a["n_facies"][]
        fa = attributes(g["facies$(i)"])
        f = Facies(
            fa["maximum_growth_rate"][] * u"m/Myr",
            fa["extinction_coefficient"][] * u"m^-1",
            fa["saturation_intensity"][] * u"W/m^2")
        depth = (0.1:0.1:50.0)u"m"
        prod = [production_rate(insolation, f, d) for d in depth]
        lines!(ax, prod / u"m/Myr", depth / u"m")
    end
end

end
# ~/~ end
