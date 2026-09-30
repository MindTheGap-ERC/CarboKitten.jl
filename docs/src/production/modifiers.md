# Production modifiers

Instead of a separate `production_modifiers` list on `Input`, time-varying
behaviour is expressed by wrapping any production spec in `MultiplyProduction`.
This implements the modifier pattern as an `AbstractProduction -> AbstractProduction`
transformation, keeping all production logic self-contained in the facies
definition.

### `MultiplyProduction`

```julia
using CarboKitten.Production: MultiplyProduction, BenthicProduction

# Reef growth halved between 0.5 and 1.0 Myr
facies = ALCAP.Facies(
    production = MultiplyProduction(
        BenthicProduction(
            maximum_growth_rate = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        0.5;
        t_range = (0.5u"Myr", 1.0u"Myr")))
```

Modifiers compose by nesting:

```julia
MultiplyProduction(
    MultiplyProduction(base_prod, 0.5; t_range=(0u"Myr", 0.3u"Myr")),
    2.0;
    t_range=(0.8u"Myr", 1.0u"Myr"))
```

Parameters:

- `base` — any `AbstractProduction` to wrap.
- `factor::Float64` — multiplicative scaling factor.
- `t_range` — `:` (always active) or a `(t_lo, t_hi)` tuple in time units.

`MultiplyProduction` delegates `is_benthic`, `is_pelagic`, and `is_interpolated`
to its `base`, so CA participation is correctly inherited.

### How it works

`production_profile(input, p::MultiplyProduction)` calls
`production_profile(input, p.base)` and wraps the result:

``` {.julia #multiply-production-profile}
function production_profile(input::AbstractInput, p::MultiplyProduction)
    base_profile = production_profile(input, p.base)
    return function(t, w)
        f = p.t_range isa Colon || (p.t_range[1] <= t <= p.t_range[2]) ? p.factor : 1.0
        return base_profile(t, w) * f
    end
end
```

Because modifiers are baked into the production closure, `uniform_production`
and `CAProduction` contain no modifier-related code — they simply call
`capped_production(profile, t, wd, dt)`.

``` {.julia #multiply-production}
# =============================================================================
# Time-window modifier — AbstractProduction transformer
# =============================================================================

const _ProdTime     = typeof(1.0u"Myr")
const _ProdTimeSpec = Union{Colon, Tuple{_ProdTime,_ProdTime}}

"""
    MultiplyProduction(base, factor; t_range=:)

Wraps `base::AbstractProduction`, multiplying its output by `factor` during
`t_range`. Outside `t_range` the base production is unchanged.

This implements the modifier pattern as `AbstractProduction -> AbstractProduction`:
modifiers compose directly in the production spec rather than in a separate
`production_modifiers` list on `Input`.
"""
@kwdef struct MultiplyProduction <: AbstractProduction
    base::AbstractProduction
    factor::Float64
    t_range::_ProdTimeSpec = (:)
end

MultiplyProduction(base, factor::Real; kwargs...) =
    MultiplyProduction(; base=base, factor=Float64(factor), kwargs...)

is_benthic(p::MultiplyProduction)      = is_benthic(p.base)
is_pelagic(p::MultiplyProduction)      = is_pelagic(p.base)
is_interpolated(p::MultiplyProduction) = is_interpolated(p.base)
```

### `ProductionBoost`

We may use the `ProductionBoost` type to have an easier interface for `MultiplyProduction`. For example, to inhibit production at one stage and boost it at a later stage, you might write something like this:

```julia
base_production = BenthicProduction(
    maximum_growth_rate = 500u"m/Myr",
    extinction_coefficient = 0.8u"m^-1",
    saturation_intensity = 60u"W/m^2")

boosted_production = base_production * 
    ProductionBoost(0.5, (0.2u"Myr", 0.4u"Myr")) * 
    ProductionBoost(1.5, (0.6u"Myr", 0.8u"Myr"))
```

The implementation is smaller than the use case:

``` {.julia #production-boost}
@kwdef struct ProductionBoost
    factor::Float64
    t_range::_ProdTimeSpec = (:)
end

Base.:*(p::AbstractProduction, b::ProductionBoost) = MultiplyProduction(p, b.factor, b.t_range)
```


``` {.julia file=src/Production/Modifiers.jl}
module Modifiers
    using Unitful

    import ..Abstract: AbstractInput, AbstractProduction, production_profile, is_benthic, is_pelagic, is_interpolated

    <<multiply-production>>
    <<multiply-production-profile>>
    <<production-boost>>
end
```
