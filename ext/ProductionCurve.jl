# ~/~ begin <<docs/src/visualization.md#ext/ProductionCurve.jl>>[init]
module ProductionCurve

using Makie
using Unitful
using HDF5

import CarboKitten.Components.Common: AbstractInput
import CarboKitten.Visualization: production_curve!, production_curve
using CarboKitten.Production: BenthicProduction, benthic_production, PelagicProduction, pelagic_production, InterpolatedProduction, production_profile
using Interpolations

function production_curve!(ax, input::I) where I <: AbstractInput
    ax.xlabel = "production [m/Myr]"
    ax.ylabel = "depth [m]"
    ax.yreversed = true

    insolation_repr = input.insolation isa Quantity ? input.insolation :
        input.insolation isa AbstractVector ? first(input.insolation) :
        input.insolation(input.time.t0)

    max_d = 50.0u"m"
    for f in input.facies
        if f.production isa InterpolatedProduction && !isempty(f.production.depth_knots)
            max_d = max(max_d, maximum(f.production.depth_knots))
        end
    end
    depth = (0.1u"m":0.1u"m":max_d)

    for (i, f) in enumerate(input.facies)
        profile = production_profile(input, f.production)
        prod = [profile(insolation_repr, d) for d in depth]
        lines!(ax, prod ./ u"m/Myr", depth ./ u"m"; label="facies $(i)")
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

    ax.xlabel = "production [m/Myr]"
    ax.ylabel = "depth [m]"

    for i in 1:a["n_facies"][]
        fa = HDF5.attributes(g["facies$(i)"])
        if fa["type"][] == "benthic"
            f = BenthicProduction(
                maximum_growth_rate = fa["maximum_growth_rate"][] * u"m/Myr",
                extinction_coefficient = fa["extinction_coefficient"][] * u"m^-1",
                saturation_intensity = fa["saturation_intensity"][] * u"W/m^2")
            depth = (0.1u"m":0.1u"m":-max_depth)
            prod = [benthic_production(insolation, f, d) for d in depth]
            lines!(ax, prod / u"m/Myr", - depth / u"m")
        end
        if fa["type"][] == "pelagic"
            f = PelagicProduction(
                maximum_growth_rate = fa["maximum_growth_rate"][] * u"1/Myr",
                extinction_coefficient = fa["extinction_coefficient"][] * u"m^-1",
                saturation_intensity = fa["saturation_intensity"][] * u"W/m^2")
            depth = (0.1u"m":0.1u"m":-max_depth)
            prod = [pelagic_production(insolation, f, d) for d in depth]
            lines!(ax, prod / u"m/Myr", - depth / u"m")
        end
        if fa["type"][] == "interpolated"
            max_rate = fa["maximum_production"][] * u"m/Myr"
            fg = g["facies$(i)"]
            knots = fg["depth_knots"][]
            mults = Float64.(fg["multipliers"][])
            order = sortperm(knots)
            itp = linear_interpolation(knots[order], mults[order], extrapolation_bc=Flat())
            plot_depth = max(-max_depth, maximum(knots) * u"m")
            depth = (0.1u"m":0.1u"m":plot_depth)
            prod = [max_rate * itp(Unitful.ustrip(u"m", d)) for d in depth]
            lines!(ax, prod / u"m/Myr", - depth / u"m")
        end
    end
end

end
# ~/~ end
