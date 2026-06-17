# ~/~ begin <<docs/src/visualization.md#ext/ProductionCurve.jl>>[init]
module ProductionCurve

using Makie
using Unitful
using HDF5
using Interpolations

import CarboKitten.Components.Common: AbstractInput
import CarboKitten.Visualization: production_curve!, production_curve

using CarboKitten.Production: production_profile
using CarboKitten.Utility: in_units_of

# ---------------------------------------------------------------------------
# In-memory form: evaluate production_profile directly
# ---------------------------------------------------------------------------

"""
    production_curve!(ax, input)

Plot production rate vs depth for each facies in `input`. The production
profile `(t, w) -> rate` is evaluated at the initial time (t=0) over a depth
range determined by the facies definitions.
"""
function production_curve!(ax, input::I) where {I <: AbstractInput}
    ax.xlabel = "production [m/Myr]"
    ax.ylabel = "depth [m]"
    ax.yreversed = true

    t0 = input.time.t0

    max_d = 50.0u"m"
    for f in input.facies
        p = f.production
        if hasproperty(p, :depth_knots) && !isempty(p.depth_knots)
            max_d = max(max_d, maximum(p.depth_knots))
        end
    end

    depth = (0.1u"m":0.1u"m":max_d)

    for (i, f) in enumerate(input.facies)
        profile = production_profile(input, f.production)
        prod = [profile(t0, d) for d in depth]
        lines!(ax, prod .|> in_units_of(u"m/Myr"), depth .|> in_units_of(u"m");
               label = "facies $(i)")
    end
end

function production_curve(input::I) where {I <: AbstractInput}
    fig = Figure()
    ax  = Axis(fig[1, 1])
    production_curve!(ax, input)
    Legend(fig[1, 2], ax)
    return fig
end

# ---------------------------------------------------------------------------
# HDF5 form: read the generic 2D production table
# ---------------------------------------------------------------------------

"""
    production_curve!(ax, g::HDF5.Group; n_time_samples=5)

Plot production rate vs depth from the 2D production table saved in the HDF5
`input` group. Each facies is shown as a family of `n_time_samples` curves
at representative times (darker = later). No knowledge of the underlying
production type is required.
"""
function production_curve!(ax, g::HDF5.Group; n_time_samples::Int = 5)
    ax.xlabel = "production [m/Myr]"
    ax.ylabel = "depth [m]"
    ax.yreversed = true

    a = HDF5.attributes(g)
    n_facies = a["n_facies"][]

    for i in 1:n_facies
        key = "facies$(i)"
        haskey(g, key) || continue
        fg = g[key]
        haskey(fg, "production_table") || continue

        table      = read(fg["production_table"])      # (n_depth, n_time)
        depth_axis = read(fg["production_depth_axis"]) # m
        n_d, n_t   = size(table)

        t_idx = unique(round.(Int, LinRange(1, n_t, min(n_time_samples, n_t))))

        for (k, ti) in enumerate(t_idx)
            alpha = 0.3 + 0.7 * k / length(t_idx)
            label = k == 1 ? "facies $(i)" : nothing
            lines!(ax, table[:, ti], depth_axis;
                   label = label, alpha = alpha)
        end
    end
end

function production_curve(filename::AbstractString; kwargs...)
    h5open(filename, "r") do fid
        fig = Figure()
        ax  = Axis(fig[1, 1])
        production_curve!(ax, fid["input"]; kwargs...)
        Legend(fig[1, 2], ax)
        return fig
    end
end

end
# ~/~ end
