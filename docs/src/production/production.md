# Production

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

## Insolation curve

Production profiles are now functions of `(time, water_depth)` rather than
`(insolation, water_depth)`.

The `insolation_curve` helper captures insolation inside the closure returned by `production_profile`, so the model loop never needs to call an insolation function explicitly.

``` {.julia #insolation-curve}
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
```

## Interface

``` {.julia file=src/Production/Abstract.jl}
module Abstract
    
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

``` {.julia file=src/Production.jl}
module Production

include("Production/Abstract.jl")
include("Production/Benthic.jl")
include("Production/Pelagic.jl")
include("Production/Interpolated.jl")
include("Production/Modifiers.jl")

using Unitful
using ..Utility: in_units_of

import .Abstract: AbstractProduction, production_profile

export AbstractProduction, production_profile

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
