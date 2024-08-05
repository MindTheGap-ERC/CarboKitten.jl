# ~/~ begin <<docs/src/denudation/empirical.md#src/Denudation/EmpiricalDenudationMod.jl>>[init]
## this code used the emperical equations with form of D = f(precipitation, slope) to estimate the denudation rates on exposed carbonate platform.
module EmpiricalDenudationMod

import ..Abstract: DenudationType, denudation, redistribution
using ...Config: Box
using Unitful

@kwdef struct EmpiricalDenudation <: DenudationType
    precip
end

#calculate planar slopes based on [ARCGIS apporach](https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/how-slope-works.htm)
function slope_kernel(w::Any, cellsize::Float64)
    dzdx = (-w[1, 1] - 2 * w[2, 1] - w[3, 1] + w[1, 3] + 2 * w[2, 3] + w[3, 3]) / (8 * cellsize)
    dzdy = (-w[1, 1] - 2 * w[1, 2] - w[1, 3] + w[3, 1] + 2 * w[3, 2] + w[1, 1]) / (8 * cellsize)

    if abs(w[2, 2]) <= abs(min(w...))
        return 0.0
    else
        atan(sqrt(dzdx^2 + dzdy^2)) * (180 / π)
    end
end
#small bug: slope not equals to 0 for depression.

# calculate denudation based on regressed function
function empirical_denudation(precip::Float64, slope::Any)
    local a = 9.1363
    local b = -0.008519
    local c = 580.51
    local d = 9.0156
    local e = -0.1245
    local f = 4.91086
    (a ./ (1 .+ exp.(b .* (precip .- c)))) .* (d ./ (1 .+ exp.(e .* (slope .- f)))) # using logistic function
end

function denudation(::Box, p::EmpiricalDenudation, water_depth, slope, facies)
    return (empirical_denudation.(p.precip, slope) .* u"m/kyr")
end

function redistribution(box::Box, ::EmpiricalDenudation, water_depth, slope, facies)
    return nothing
end

end
# ~/~ end
