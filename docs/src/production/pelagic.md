# Pelagic Production

A facies can be specified as being pelagic, meaning that production is not governed at the sea floor (benthic zone), rather the entire water column above contributes to the production. The local production rate still follows the same function as before (well motivated by exponential decay of available radiation), but now we need to integrate over the water depth:

$$g_p(w) = \int_0^{w} g_m \tanh\left({{I_0 e^{-kw}} \over {I_k}}\right) \textrm{d}w.$$

This integral has no analytic solution, so we'll use a numeric integrator to evaluate the complete production curve, and then linearly interpolate the generated table to obtain production rates.

Pelagic facies do not participate in the CA.

``` {.julia #pelagic-production}
function pelagic_production(insolation, facies, water_depth)
    return quadgk(w -> production_rate(insolation, facies, w), 0.0u"m", water_depth)[1]
end
```

``` {.julia file=examples/production/pelagic.jl}
module PelagicProductionPlot

using CarboKitten
using CarboKitten.Production
using CairoMakie

@kwdef struct Input <: CarboKitten.AbstractInput
    insolation = 400.0u"W/m^2"
end

function main()
    water_depth = (0.01:0.1:50.0)u"m"
    fig = Figure(size=(600,600))
    input = Input()

    ax = Axis(fig[1, 1], yreversed=true, ylabel="depth [m]", xlabel="production [m/Myr]")
    for (k, prod) in pairs(Production.EXAMPLE)
        f = production_profile(input, prod)
        p = water_depth .|> (w -> f(input.time.t0, w))
        lines!(ax, p |> in_units_of(u"m/Myr"),
            water_depth |> in_units_of(u"m"), label = string(k))
    end
    fig[1, 2] = Legend(fig, ax, "Facies")
    fig
end

end
```
