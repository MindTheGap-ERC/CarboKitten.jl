# ~/~ begin <<docs/src/visualization.md#ext/ProductionCurve.jl>>[init]
module ProductionCurve

using Makie
using Unitful
using HDF5

import CarboKitten.Components.Common: AbstractInput
import CarboKitten.Visualization: production_curve!, production_curve
using CarboKitten.Components.Production: Facies, production_rate

function production_curve!(ax, input::I) where I <: AbstractInput
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

function production_curve(input::I) where I <: AbstractInput
    fig = Figure()
    ax = Axis(fig[1, 1])
    production_curve!(ax, input)
    fig
end

function production_curve(filename::AbstractString)
    h5open(filename, "r") do fid
        fig = Figure()
        ax = Axis(fig[1, 1])
        production_curve!(ax, fid["input"])
        fig
    end
end

function production_curve!(ax, g::HDF5.Group; max_depth=-50.0u"m")
    a = HDF5.attributes(g)
    insolation = 400.0u"W/m^2"  # a["insolation"][] * u"W/m^2"

    ax.title = "production at $(sprint(show, insolation; context=:fancy_exponent=>true))"
    ax.xlabel = "production [m/Myr]"
    ax.ylabel = "depth [m]"

    for i in 1:a["n_facies"][]
        fa = HDF5.attributes(g["facies$(i)"])
        f = Facies(
            maximum_growth_rate = fa["maximum_growth_rate"][] * u"m/Myr",
            extinction_coefficient = fa["extinction_coefficient"][] * u"m^-1",
            saturation_intensity = fa["saturation_intensity"][] * u"W/m^2")
        depth = (0.1u"m":0.1u"m":-max_depth)
        prod = [production_rate(insolation, f, d) for d in depth]
        lines!(ax, prod / u"m/Myr", - depth / u"m")
    end
end

end
# ~/~ end