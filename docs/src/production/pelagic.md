# Pelagic Production

A facies can be specified as being pelagic, meaning that production is not governed at the sea floor (benthic zone), rather the entire water column above contributes to the production. The local production rate still follows the same function as before (well motivated by exponential decay of available radiation), but now we need to integrate over the water depth:

```math
_p(w) = \int_0^{w} g_m \tanh\left({{I_0 e^{-kw}} \over {I_k}}\right) \textrm{d}w.
```

This integral has no analytic solution, so we'll use a numeric integrator to evaluate the complete production curve, and then linearly interpolate the generated table to obtain production rates.

Pelagic facies do not participate in the CA.

``` {.julia #pelagic-production}
function pelagic_production(insolation, facies, water_depth)
    return quadgk(w -> production_rate(insolation, facies, w), 0.0u"m", water_depth)[1]
end
```

Because the parameters for benthic and pelagic production have different units, we need different types to store them.

``` {.julia file=src/Production/Pelagic.jl}
module Pelagic

using Unitful
using Interpolations
using QuadGK

using ...Utility: in_units_of
import ..Abstract: AbstractInput, AbstractProduction, is_pelagic, production_profile
import ..Insolation: insolation_curve
import ..Benthic: production_rate

@kwdef struct PelagicProduction <: AbstractProduction
    maximum_growth_rate::typeof(1.0u"1/Myr") = 0.0u"1/Myr"
    extinction_coefficient::typeof(1.0u"m^-1") = 0.0u"m^-1"
    saturation_intensity::typeof(1.0u"W/m^2") = 1.0u"W/m^2"
    maximum_production_depth::typeof(1.0u"m") = 200.0u"m"
    table_size::Tuple{Int, Int} = (1000, 1000)
end

is_pelagic(::PelagicProduction) = true

<<pelagic-production>>
<<production-lookup>>

production_profile(input::AbstractInput, p::PelagicProduction) =
    pelagic_production_lookup(input, p)

end
```

## Lookup tables

We use `linear_interpolation` from `Interpolations` to compute production profiles from look-up tables. `insolation_curve` provides a `time -> insolation` closure used to evaluate the lookup at the correct insolation for each time step.

``` {.julia #production-lookup}
function pelagic_production_lookup(input::AbstractInput, prod::PelagicProduction)
    I_of_t = insolation_curve(input)
    depth_grid = LinRange(0.0, prod.maximum_production_depth |> in_units_of(u"m"), prod.table_size[2])

    if input.insolation isa Quantity
        # Constant insolation — 1D depth lookup, time argument ignored
        I0 = input.insolation
        production_values = [pelagic_production(I0, prod, w * u"m") |> in_units_of(u"m/Myr")
                             for w in depth_grid]
        itp = linear_interpolation(depth_grid, production_values, extrapolation_bc=Line())
        return (_, w) -> itp(w |> in_units_of(u"m")) * u"m/Myr"
    end

    # Variable insolation — 2D (insolation × depth) lookup
    t_axis = time_axis(input)
    I_vals = [I_of_t(t) |> in_units_of(u"W/m^2") for t in t_axis]
    I_min, I_max = extrema(I_vals)

    insolation_grid = LinRange(I_min, I_max, prod.table_size[1])
    production_values = [
        pelagic_production(I * u"W/m^2", prod, w * u"m") |> in_units_of(u"m/Myr")
        for I in insolation_grid, w in depth_grid
    ]
    itp = linear_interpolation((collect(insolation_grid), collect(depth_grid)),
                               production_values, extrapolation_bc=Line())
    return (t, w) -> itp(I_of_t(t) |> in_units_of(u"W/m^2"), w |> in_units_of(u"m")) * u"m/Myr"
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
