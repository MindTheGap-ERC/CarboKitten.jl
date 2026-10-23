# ~/~ begin <<docs/src/production/interpolated.md#src/Production/Interpolated.jl>>[init]
module Interpolated

using Unitful
using ...Utility: in_units_of

using Interpolations
import ..Abstract: AbstractInput, AbstractProduction, production_profile
# =============================================================================
# Interpolated (knot-based) production curve
# =============================================================================

"""
    InterpolatedProduction(; maximum_production, depth_knots, multipliers)

A depth-only production curve defined by a peak rate and a piecewise-linear
shape over `(depth, multiplier)` knots. Independent of insolation.

    rate(t, w) = maximum_production × interpolate(depth_knots, multipliers; w)

`depth_knots` need not be sorted; they are sorted internally.

# Example

    InterpolatedProduction(
        maximum_production = 500.0u"m/Myr",
        depth_knots        = [0.0u"m", 5.0u"m", 20.0u"m", 50.0u"m"],
        multipliers        = [0.0,     1.0,     0.6,      0.0])
"""
@kwdef struct InterpolatedProduction <: AbstractProduction
    maximum_production::typeof(1.0u"m/Myr") = 0.0u"m/Myr"
    depth_knots::Vector{typeof(1.0u"m")}    = typeof(1.0u"m")[]
    multipliers::Vector{Float64}            = Float64[]
end

is_benthic(::InterpolatedProduction)      = false
is_pelagic(::InterpolatedProduction)      = false
is_interpolated(::InterpolatedProduction) = true

function production_profile(::AbstractInput, p::InterpolatedProduction)
    @assert length(p.depth_knots) == length(p.multipliers)
    @assert length(p.depth_knots) >= 2
    depths_m = [d |> in_units_of(u"m") for d in p.depth_knots]
    order = sortperm(depths_m)
    itp = linear_interpolation(depths_m[order], p.multipliers[order], extrapolation_bc=Flat())
    max_rate = p.maximum_production
    return (_, w) -> max_rate * itp(w |> in_units_of(u"m"))
end

end
# ~/~ end
