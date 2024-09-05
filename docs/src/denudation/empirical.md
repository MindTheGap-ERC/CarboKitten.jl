# Emperical denudation
Cl isotopes are an emerging tool to decipher the denudation rates (chemical dissolution + physical erosion) in carbonate-dominated area.

Research based on the karst region and carbonate platform terrace suggested that the denudation rates are mainly controlled by precipitation and slopes, although the debates about which factor is more important is still ongoing ([yang_combined_2020](@cite), [thomas_limited_2018](@cite)). In general, the precipitation mainly controls the chemical dissolution while the slopes mainly controls the physical ersions. In addition, the type of carbonates may also play an important role ([krklec_long-term_2022](@cite)), but given this feature is studied poorly so we will ditch it for now. We have checked and compiled the denudation rates (mm/kyr) along with precipitation and slopes serve as a starting point to create a function relates denudation rates (mm/kyr) to precipitation and slopes. The compiled data could be found in OSFdatabase. This is an empirical relationship and have a relatively large uncertainty in terms of fitting.

<img width="433" alt="image" src="https://github.com/MindTheGap-ERC/CarboKitten.jl/assets/64159957/29da37f3-0e14-479f-8ed3-a8b5adc052a5">

*Fig 1. The relationship between MAP (mean precipitation per year, mm/y) and denudation rates (mm/ky)*

<img width="370" alt="image" src="https://github.com/MindTheGap-ERC/CarboKitten.jl/assets/64159957/1bd9e2ae-d4d7-4611-a8b8-83459222022f">

*Fig 2. The relationship between the curvature (i.e., slope or height) and the denudation rates (mm/ky)*

We can see that both the slope and precipitation could increase the denudation rates, and reaches a 'steady state' after a certain point.

In terms of implementation, I used the function form of $D = P * S$, where $D$ means denudation rates, $P$ means effects of precipitation while $S$ means effects of Slope. By doing so, we can consider both effects. Such formula structure is similar to RUSLE model, a widely used LEM (e.g., (thapa_spatial_2020)[@cite]).

``` {.julia file=src/Denudation/EmpiricalDenudationMod.jl}
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
    (a ./ (1 .+ exp.(b .* (precip .* 1000 .- c)))) .* (d ./ (1 .+ exp.(e .* (slope .- f)))) .* u"m/kyr" # using logistic function
end

function denudation(::Box, p::EmpiricalDenudation, water_depth, slope, facies)
    precip = p.precip ./ u"m"
    return empirical_denudation.(precip, slope)
end

function redistribution(::Box, p::EmpiricalDenudation, denudation_mass, water_depth)
    return denudation_mass .* 0
end

end
```
