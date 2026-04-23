# ~/~ begin <<docs/src/visualization.md#ext/ProductionCurve.jl>>[init]
module ProductionCurve

using Makie
using Unitful
using HDF5

import CarboKitten.Components.Common: AbstractInput
import CarboKitten.Visualization: production_curve!, production_curve, facies_colormap
using CarboKitten.Components.Production: Facies, production_rate



function production_curve!(
    ax,
    input::I;
    facies_colors = nothing,
    max_depth = 500.0u"m",
    facies_names = nothing
) where I <: AbstractInput

    ax.title = "production at $(sprint(show, input.insolation; context=:fancy_exponent=>true))"
    ax.xlabel = "production [m/Myr]"
    ax.ylabel = "depth [m]"
    ax.yreversed = true

    t0 = input.time.t0
    n_f = length(input.facies)

    # fallback names
    if facies_names === nothing
        facies_names = string.(1:n_f)
    end

    cmap = facies_colormap(n_f; facies_colors=facies_colors, include_nodeposit=false)
    colors = Makie.to_color.(cmap.colors)

    depth = (0.1u"m":0.5u"m":max_depth)

    for (i, f) in enumerate(input.facies)
        prod = [production_rate(f, d, t0) for d in depth]

        lines!(ax,
            prod / u"m/Myr",
            depth / u"m";
            color = colors[i],
            label = facies_names[i]
        )
    end

    axislegend(ax)
end


function production_curve(
    input::I;
    facies_colors = nothing,
    max_depth = 500.0u"m",
    facies_names = nothing
) where I <: AbstractInput

    fig = Figure()
    ax = Axis(fig[1,1])

    production_curve!(
        ax,
        input;
        facies_colors = facies_colors,
        max_depth = max_depth,
        facies_names = facies_names
    )

    return fig
end


function production_curve!(ax, g::HDF5.Group; max_depth=-50.0u"m", facies_colors=nothing, facies_names = nothing)
    a = HDF5.attributes(g)
    insolation = 400.0u"W/m^2"  # a["insolation"][] * u"W/m^2"

    ax.title = "production at $(sprint(show, insolation; context=:fancy_exponent=>true))"
    ax.xlabel = "production [m/Myr]"
    ax.ylabel = "depth [m]"
	ax.yreversed = true
    n_f = a["n_facies"][]
cmap = facies_colormap(n_f; facies_colors=facies_colors, include_nodeposit=false)
colors = Makie.to_color.(cmap.colors)

for i in 1:n_f
        fa = HDF5.attributes(g["facies$(i)"])
        gmax = fa["maximum_growth_rate"][] * u"m/Myr"

depth_knots = Tuple{typeof(1.0u"m"), Float64}[]
if haskey(fa, "depth_knots/depth_m")
    depths = fa["depth_knots/depth_m"][] * u"m"
    mults  = fa["depth_knots/mult"][]
    depth_knots = [(depths[j], Float64(mults[j])) for j in eachindex(mults)]
end

time_windows = Tuple{typeof(1.0u"Myr"), typeof(1.0u"Myr"), Float64}[]
if haskey(fa, "time_windows/t1_Myr")
    t1 = fa["time_windows/t1_Myr"][] * u"Myr"
    t2 = fa["time_windows/t2_Myr"][] * u"Myr"
    mm = fa["time_windows/mult"][]
    time_windows = [(t1[j], t2[j], Float64(mm[j])) for j in eachindex(mm)]
end

f = Facies(maximum_growth_rate=gmax, depth_knots=depth_knots, time_windows=time_windows)

t0 = haskey(a, "t0") ? (a["t0"][] * u"Myr") : 0.0u"Myr"
depth = (0.1u"m":0.5u"m":max_depth)
prod = [production_rate(f, d, t0) for d in depth]

        lines!(ax,
    prod / u"m/Myr",
    depth / u"m";
    color = colors[i],
    label = facies_names[i]
)

axislegend(ax)
    end
end

end
# ~/~ end