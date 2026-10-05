# Production

CarboKitten supports a variety of methods to model sediment production. In general a production model is a function of waterdepth and insolation giving a value in units of meters (of sediment) per million years. We provide two models based on Bosscher & Schlager 1992 [Bosscher1992](@cite): `BenthicProduction` and `PelagicProduction`.

## Benthic and Pelagic production models

The *benthic* production model supposes that organisms (like corals or bivalves) on the bottom of the sea are responsible for producing sediment, depending on how much sun light is able to reach the bottom. At a given saturation intensity, the organism is not able to produce any more sediment, due to other limiting factors. The growth rate is given as,

```math
g_b(w) = g_m \tanh\left({{I_0 e^{-kw}} \over {I_k}}\right).
```

In the case of *pelagic* production, the supposed organism (maybe algae) occupy the entire water column above the sea floor. In that case we need to integrate the benthic production profile over the water column to obtain the pelagic production profile,

```math
g_p(w) = \int_0^{w} g_m \tanh\left({{I_0 e^{-kw}} \over {I_k}}\right) \textrm{d}w.
```

Both these models have the same parameters, albeit that ``g_m`` has slightly different units:

- ``g_m`` is the maximum production rate in units of ``{\rm m/Myr}`` (or ``{\rm 1/Myr}`` in case of pelagic production).
- ``I_0`` is the given insolation in ``{\rm W}/{\rm m}^2`` (set as a separate global input paramater, possibly a function of time).
- ``I_k`` is the saturation intensity in ``{\rm W}/{\rm m}^2``.
- ``k`` is the extinction coefficient in ``1/{\rm m}``.
- ``w`` is the water depth as obtained from CarboKitten's model state.

The implementation of these models can be found in their respective sections: [Benthic Production](@ref) and [Pelagic Production](@ref).

## Interpolated and modified production

As stated before, a production model is a function of waterdepth and insolation giving a value in units of meters (of sediment) per million years. More abstractly, we can replace the insolation with a time dependency to allow more flexible time dependent production levels. We should warn that, in general, more flexibility does impact the predictive power of the model negatively.

The [Interpolated Production](@ref) model lets the user provide a set of key values for production as a function of water depth, while a separate set of [Production Modifiers](@ref) can be used to modulate production amplitudes over time. These modifiers can be chained and applied to any other production model.

## No production

In some cases you may whish to disable production altogether. For that case we have defined the `NoProduction` model. This is the default setting for the `Facies` configuration.

## Mixing curve types

Different facies in the same run can use different production types:

```julia
facies = [
    ALCAP.Facies(production = BenthicProduction(
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2")),
    ALCAP.Facies(production = InterpolatedProduction(
        maximum_production=350.0u"m/Myr",
        depth_knots=[0.0u"m", 10.0u"m", 40.0u"m", 60.0u"m"],
        multipliers=[0.0, 0.5, 1.0, 0.0])),
]
```

## Custom production models

If none of the above options provide what you need, you may choose to implement your own production model. Derive a new `struct` from `AbstractProduction`, and implement the `production_profile` method:

```julia
module MyProductionModel
    
import CarboKitten: production_profile, AbstractProduction, AbstractInput, insolation_curve

struct MyProduction <: AbstractProduction
    ...
end

production_profile(input::AbstractInput, production::MyProduction) =
    I_of_t = insolation_curve(input)
    function (time, water_depth)
        I = I_of_t(time)
        ...
    end
    
end  # module MyProductionModel
```

### Insolation curve

The `insolation_curve` helper captures insolation inside the closure returned by `production_profile`, so the model loop never needs to call an insolation function explicitly.

``` {.julia file=src/Production/Insolation.jl}
module Insolation

using Unitful

import ..Abstract: AbstractInput

"""
    insolation_curve(input) -> function(time) -> insolation

Build a closure mapping time to insolation from the input specification.
Handles constant (`Quantity`), tabular (`AbstractVector`), and functional
insolation inputs.
"""
function insolation_curve(input::AbstractInput)
    insolation_param = input.insolation
    if insolation_param isa Quantity
        return _ -> insolation_param
    elseif insolation_param isa AbstractVector
        t_axis = time_axis(input)
        t_vals = ustrip.(u"Myr", t_axis)
        I_vals = ustrip.(u"W/m^2", insolation_param)
        itp = linear_interpolation(t_vals, I_vals, extrapolation_bc=Flat())
        return t -> itp(ustrip(u"Myr", t)) * u"W/m^2"
    else
        return t -> insolation_param(t)
    end
end

end
```

## Interface

Every production model should adhere to the production interface outlined below.

``` {.julia file=src/Production/Abstract.jl}
module Abstract

using Unitful
using ...CarboKitten: AbstractInput

"""
    production_profile(input::AbstractInput, p)

Given an input and a production configuration, returns a function
`(time, water_depth) -> production_rate`.

Insolation is read from `input` and composed into the returned closure via
`insolation_curve` — the caller does not need to pass insolation at each
time step.

The default implementation assumes `p` is already a callable
`(time, water_depth) -> rate`.

## Example

    struct MyProduction <: AbstractProduction
        ...
    end

    import CarboKitten.Production: production_profile

    production_profile(input::AbstractInput, p::MyProduction) =
        function (time, water_depth)
            ...
        end
"""
production_profile(input::AbstractInput, p) = p

"""
    is_benthic(obj)

Predicate to determine if a facies or production spec is benthic.
Defaults to `false`.
"""
is_benthic(p) = false

"""
    is_pelagic(obj)

Predicate to determine if a facies or production spec is pelagic.
Defaults to `false`.
"""
is_pelagic(p) = false

"""
    is_interpolated(obj)

Predicate to determine if a facies or production spec is interpolation-based.
Defaults to `false`.
"""
is_interpolated(p) = false

abstract type AbstractProduction end

struct NoProduction <: AbstractProduction
end

production_profile(::AbstractInput, ::NoProduction) = (_, _) -> 0.0u"m/Myr"

end
```

## Module

The production module collects the different production models, and provides a list of example profiles.

``` {.julia file=src/Production.jl}
module Production

include("Production/Abstract.jl")
include("Production/Insolation.jl")
include("Production/Benthic.jl")
include("Production/Pelagic.jl")
include("Production/Interpolated.jl")
include("Production/Modifiers.jl")

using Unitful

import .Abstract: AbstractProduction, NoProduction, production_profile, is_benthic, is_pelagic, is_interpolated
import .Benthic: BenthicProduction
import .Pelagic: PelagicProduction
import .Interpolated: InterpolatedProduction
import .Modifiers: MultiplyProduction, ProductionBoost
import .Insolation: insolation_curve

export AbstractProduction, production_profile, NoProduction, BenthicProduction, PelagicProduction,
    InterpolatedProduction, MultiplyProduction, ProductionBoost, EXAMPLE

const EXAMPLE = Dict(
    :euphotic => BenthicProduction(
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2"),
    :oligophotic => BenthicProduction(
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2"),
    :aphotic => BenthicProduction(
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2"),
    :pelagic => PelagicProduction(
        maximum_growth_rate=7.0u"1/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2"),
    :interpolated => InterpolatedProduction(
        maximum_production=500u"m/Myr",
        depth_knots=[0.0u"m", 5.0u"m", 15.0u"m", 30.0u"m", 50.0u"m"],
        multipliers=[0.0, 1.0, 1.0, 0.4, 0.0]),
    :time_varying => MultiplyProduction(
        BenthicProduction(
            maximum_growth_rate=500u"m/Myr",
            extinction_coefficient=0.8u"m^-1",
            saturation_intensity=60u"W/m^2"),
        0.5;
        t_range=(0.5u"Myr", 1.0u"Myr"))
)

end
```
