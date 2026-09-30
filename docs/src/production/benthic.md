# Benthic Production

The `Production` module specifies the production rate following the model by Bosscher & Schlager 1992 [Bosscher1992](@cite).
The growth rate is given as

```math
g(w) = g_m \tanh\left({{I_0 e^{-kw}} \over {I_k}}\right).
```

This can be understood as a smooth transition between the maximum growth rate under saturated conditions, and exponential decay due to light intensity dropping with greater water depth.

``` {.julia #benthic-production-rate}
function production_rate(insolation, facies, water_depth)
    gₘ = facies.maximum_growth_rate
    I = insolation / facies.saturation_intensity
    x = water_depth * facies.extinction_coefficient
    return x > 0.0 ? gₘ * tanh(I * exp(-x)) : zero(typeof(gₘ))
end
```

We also have an alias for benthic production rates `benthic_production`, as we will also allow for a `pelagic_production` function (see [Pelagic Production](@ref)).

``` {.julia #benthic-production-rate}
benthic_production(i, f, w) = production_rate(i, f, w)
```

Insolation is captured inside each production profile closure via `insolation_curve` — the model loop only needs to pass the current simulation time.

From just this equation we can define a uniform production process. This requires that we have a `Facies` that defines the `maximum_growth_rate`, `extinction_coefficient` and `saturation_intensity`.

The `insolation` input may be given as a scalar quantity, say `400u"W/m^2"`, or as a function of time.

``` {.julia file=src/Production/Benthic.jl}
module Benthic

using Unitful
import ..Abstract: AbstractInput, AbstractProduction, is_benthic, production_profile
import ..Insolation: insolation_curve

<<benthic-production-rate>>

@kwdef struct BenthicProduction <: AbstractProduction
    maximum_growth_rate::typeof(1.0u"m/Myr") = 0.0u"m/Myr"
    extinction_coefficient::typeof(1.0u"m^-1") = 0.0u"m^-1"
    saturation_intensity::typeof(1.0u"W/m^2") = 1.0u"W/m^2"
end

is_benthic(::BenthicProduction) = true

function production_profile(input::AbstractInput, p::BenthicProduction)
    I_of_t = insolation_curve(input)
    return (t, w) -> benthic_production(I_of_t(t), p, w)
end

end
```
