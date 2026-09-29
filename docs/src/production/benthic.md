# Benthic Production

The `Production` module specifies the production rate following the model by Bosscher & Schlager 1992 [Bosscher1992](@cite).
The growth rate is given as

$$g(w) = g_m \tanh\left({{I_0 e^{-kw}} \over {I_k}}\right).$$

This can be understood as a smooth transition between the maximum growth rate under saturated conditions, and exponential decay due to light intensity dropping with greater water depth.

``` {.julia #component-production-rate}
function production_rate(insolation, facies, water_depth)
    gₘ = facies.maximum_growth_rate
    I = insolation / facies.saturation_intensity
    x = water_depth * facies.extinction_coefficient
    return x > 0.0 ? gₘ * tanh(I * exp(-x)) : zero(typeof(gₘ))
end
```

We also have an alias for benthic production rates `benthic_production`, as we will also allow for a `pelagic_production` function (see [Pelagic Production](@ref)).

``` {.julia #component-production-rate}
benthic_production(i, f, w) = production_rate(i, f, w)
```

Because we can only produce as much as keeps the factory submerged, we have to cap the total production in a single time step to the current water depth. This assumes we have production as a function of time and water depth.

``` {.julia #component-production-rate}
"""
    capped_production(f, time, water_depth, dt)

Apply production function `f(time, water_depth) -> rate`, clip to non-negative,
and cap by available accommodation. Returns the deposited thickness for `dt`.
"""
function capped_production(f, time, water_depth, dt)
    clip_positive(x::T) where {T} = max(x, zero(T))
    p = clip_positive(f(time, water_depth))
    return min(max(0.0u"m", water_depth), p * dt)
end
```

Insolation is captured inside each production profile closure via `insolation_curve` — the model loop only needs to pass the current simulation time.

From just this equation we can define a uniform production process. This requires that we have a `Facies` that defines the `maximum_growth_rate`, `extinction_coefficient` and `saturation_intensity`.

The `insolation` input may be given as a scalar quantity, say `400u"W/m^2"`, or as a function of time.
