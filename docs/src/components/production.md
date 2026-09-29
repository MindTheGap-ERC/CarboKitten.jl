# Production

```component-dag
CarboKitten.Components.Production
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

### Mixing curve types

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

## Production Component

### Production Properties

Because the parameters for benthic and pelagic production have different units, we need different types to store them.

``` {.julia #production-profile}
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

@kwdef struct BenthicProduction <: AbstractProduction
    maximum_growth_rate::typeof(1.0u"m/Myr") = 0.0u"m/Myr"
    extinction_coefficient::typeof(1.0u"m^-1") = 0.0u"m^-1"
    saturation_intensity::typeof(1.0u"W/m^2") = 1.0u"W/m^2"
end

is_benthic(::BenthicProduction) = true
is_pelagic(::BenthicProduction) = false

function production_profile(input::AbstractInput, p::BenthicProduction)
    I_of_t = insolation_curve(input)
    return (t, w) -> benthic_production(I_of_t(t), p, w)
end

@kwdef struct PelagicProduction <: AbstractProduction
    maximum_growth_rate::typeof(1.0u"1/Myr") = 0.0u"1/Myr"
    extinction_coefficient::typeof(1.0u"m^-1") = 0.0u"m^-1"
    saturation_intensity::typeof(1.0u"W/m^2") = 1.0u"W/m^2"
    maximum_production_depth::typeof(1.0u"m") = 200.0u"m"
    table_size::Tuple{Int, Int} = (1000, 1000)
end

is_benthic(::PelagicProduction) = false
is_pelagic(::PelagicProduction) = true
production_profile(input::AbstractInput, p::PelagicProduction) = pelagic_production_lookup(input, p)
```

### Input parameters

``` {.julia #production-input}
@kwdef struct Input <: AbstractInput
    insolation
end

@kwdef struct Facies <: AbstractFacies
    production = NoProduction()
end

is_benthic(facies::AbstractFacies) = is_benthic(facies.production)
is_pelagic(facies::AbstractFacies) = is_pelagic(facies.production)
```

The `production_modifiers` field has been removed from `Input`. Time-dependent
behaviour is now expressed by composing production specs directly in the `Facies`
definition using `MultiplyProduction`.

### Pelagic Lookup tables

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

### HDF5 serialization

Since `production_profile` is now fully generic, production is saved as a
2D evaluated table per facies rather than as type-specific parameters. This
means no type-specific branches are needed in either the writer or the reader.

- `input/facies_N/production_table` — 2D array of shape `(n_depth, n_time)` in `m/Myr`.
- `input/facies_N/production_depth_axis` — depth values in meters.
- `input/facies_N/production_time_axis` — time values in Myr.

The `ProductionCurve` visualisation reads the table directly and plots
production vs depth as a family of time-slice curves.

### Boilerplate

We put the basic production equations in a separate module `CarboKitten.Production`. This will contain the dynamic API (`production_profile`) and the implementation of `BenthicProduction` and `PelagicProduction`. We isolate these from the component `CarboKitten.Components.Production` implementation to prevent these production objects from being replicated in derived components, leading to problems with dispatch on `production_profile`.

``` {.julia file=src/Production.jl}
module Production

using Unitful
using QuadGK
using Interpolations

import ..CarboKitten: AbstractInput, time_axis
using ..Utility: in_units_of

export production_profile

<<component-production-rate>>
<<pelagic-production>>
<<insolation-curve>>
<<production-profile>>
<<production-lookup>>
<<interpolated-production>>
<<multiply-production>>

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

``` {.julia file=src/Components/Production.jl}
@compose module Production
@mixin TimeIntegration, WaterDepth, FaciesBase
using ..Common
using ..WaterDepth: water_depth
using ..TimeIntegration: time, write_times
using ...Production: NoProduction, InterpolatedProduction, MultiplyProduction
import ...Production: production_profile, is_benthic, is_pelagic, is_interpolated,
    capped_production, insolation_curve

using HDF5
using QuadGK
using Interpolations
using Logging

export uniform_production
export MultiplyProduction, InterpolatedProduction, NoProduction

<<production-input>>

function write_header(input::AbstractInput, output::AbstractOutput)
    # Insolation time series
    I_of_t = insolation_curve(input)
    t_write = write_times(input)[1:end-1]
    set_attribute(output, "insolation",
        [I_of_t(t) |> in_units_of(u"W/m^2") for t in t_write])

    # Generic 2D production table — no type-specific branches needed
    depth_grid = LinRange(0.0u"m", 200.0u"m", 200)
    for (i, f) in enumerate(input.facies)
        prof = production_profile(input, f.production)
        table = [prof(t, d) |> in_units_of(u"m/Myr")
                 for d in depth_grid, t in t_write]
        set_attribute(output, "facies$(i)/production_table", table)
        set_attribute(output, "facies$(i)/production_depth_axis",
            collect(depth_grid) .|> in_units_of(u"m"))
        set_attribute(output, "facies$(i)/production_time_axis",
            t_write .|> in_units_of(u"Myr"))
    end
end

function uniform_production(input::AbstractInput)
    w = water_depth(input)
    na = [CartesianIndex()]
    facies = input.facies
    dt = input.time.Δt
    production_rates = [production_profile(input, f.production) for f in facies]
    get_time = time(input)

    p(state::AbstractState, wd::AbstractMatrix) = begin
        t = get_time(state)
        capped_production.(production_rates[:, na, na], t, wd[na, :, :], dt)
    end
    p(state::AbstractState) = p(state, w(state))
    return p
end

end
```

## CA Production

```component-dag
CarboKitten.Components.CAProduction
```

The `CAProduction` component gives production that depends on the provided CA.
Insolation is captured inside each production spec closure — the production
loop only needs the current time, not an insolation value.

``` {.julia file=src/Components/CAProduction.jl}
@compose module CAProduction
    @mixin TimeIntegration, CellularAutomaton, Production
    using ..Common
    using ..TimeIntegration: time
    using ..WaterDepth: water_depth
    using ...Production: production_profile, capped_production
    using Logging

    function production(input::AbstractInput)
        w = water_depth(input)
        na = [CartesianIndex()]
        output_ = Array{Amount, 3}(undef, n_facies(input), input.box.grid_size...)

        facies = input.facies
        dt = input.time.Δt
        production_specs = ((production_profile(input, f.production) for f in facies)...,)
        get_time = time(input)

        function p(state::AbstractState, wd::AbstractMatrix)::Array{Amount,3}
            output::Array{Amount, 3} = output_
            t = get_time(state)
            for i in eachindex(IndexCartesian(), wd)
                for f in eachindex(facies)
                    if facies[f].active
                        output[f, i[1], i[2]] = f != state.ca[i] ? 0.0u"m" :
                            capped_production(production_specs[f], t, wd[i], dt)
                    else
                        output[f, i[1], i[2]] =
                            capped_production(production_specs[f], t, wd[i], dt)
                    end
                end
            end
            return output
        end

        @inline p(state::AbstractState) = p(state, w(state))
        return p
    end
end
```

## Tests

### Production higher in shallower water

And reversed with pelagic production. The interpolated test uses a sloping
topography so different cells have different water depths.

```{.julia #production-spec}
@testset "Components/Production" begin
    let prod = BenthicProduction(
            maximum_growth_rate = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        input = Input(
            box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
            time = TimeProperties(Δt=1.0u"kyr", steps=10),
            sea_level = t -> 0.0u"m",
            initial_topography = (x, y) -> -10u"m",
            subsidence_rate = 0.0u"m/Myr",
            facies = [Facies(production=prod)],
            insolation = 400.0u"W/m^2")

        state = initial_state(input)
        prod = uniform_production(input)(state)
        @test all(prod[1:end-1,:] .>= prod[2:end,:])
    end

    let prod = PelagicProduction(
            maximum_growth_rate = 5u"1/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        input = Input(
            box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
            time = TimeProperties(Δt=1.0u"kyr", steps=10),
            sea_level = t -> 0.0u"m",
            initial_topography = (x, y) -> -10u"m",
            subsidence_rate = 0.0u"m/Myr",
            facies = [Facies(production=prod)],
            insolation = 400.0u"W/m^2")

        state = initial_state(input)
        prod = uniform_production(input)(state)
        @test all(prod[1:end-1,:] .<= prod[2:end,:])
    end
end

@testset "Components/Production/interpolated" begin
    let prod = InterpolatedProduction(
            maximum_production = 500u"m/Myr",
            depth_knots = [0.0u"m", 5.0u"m", 15.0u"m", 50.0u"m"],
            multipliers = [0.0,     1.0,     1.0,      0.0]),
        input = Input(
            box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=5.0u"m"),
            time = TimeProperties(Δt=1.0u"kyr", steps=10),
            sea_level = t -> 0.0u"m",
            initial_topography = (x, y) -> -x * 1.0,
            subsidence_rate = 0.0u"m/Myr",
            facies = [Facies(production=prod)],
            insolation = 400.0u"W/m^2")

        state = initial_state(input)
        p = uniform_production(input)(state)
        @test p[1, 3, 1] > p[1, 1, 1]
        @test all(p .>= 0.0u"m")
    end
end

@testset "Components/Production/time_varying" begin
    let base = BenthicProduction(
            maximum_growth_rate = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        prod = MultiplyProduction(base, 0.5; t_range=(0.0u"Myr", 0.5u"Myr")),
        input = Input(
            box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
            time = TimeProperties(Δt=0.1u"Myr", steps=10),
            sea_level = t -> 0.0u"m",
            initial_topography = (x, y) -> -10u"m",
            subsidence_rate = 0.0u"m/Myr",
            facies = [Facies(production=prod)],
            insolation = 400.0u"W/m^2")

        state = initial_state(input)
        # Inside t_range: production halved relative to base
        state.step = 1
        p_inside = copy(uniform_production(input)(state))
        # Outside t_range: production at full rate
        state.step = 8
        p_outside = copy(uniform_production(input)(state))
        @test all(p_outside .>= p_inside)
    end
end
```

### Variable insolation

If insolation increases linearly with time, production at t = 10 should be higher than at t = 1.

```{.julia #production-spec}
@testset "Components/Production/variable_insolation" begin
    let prod = BenthicProduction(
            maximum_growth_rate = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        input = Input(
            box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
            time = TimeProperties(Δt=1.0u"kyr", steps=10),
            sea_level = t -> 0.0u"m",
            initial_topography = (x, y) -> -10u"m",
            subsidence_rate = 0.0u"m/Myr",
            facies = [Facies(production=prod)],
            insolation = t -> 40.0u"W/m^2/kyr" * t)

        state = initial_state(input)
        state.step = 1
        prod1 = copy(uniform_production(input)(state))
        state.step = 10
        prod2 = copy(uniform_production(input)(state))
        @test all(prod2 .> prod1)
    end
end
```

``` {.julia file=test/Components/ProductionSpec.jl}
module ProductionSpec
    using Test
    using CarboKitten
    using CarboKitten.Components.Common
    using CarboKitten.Components.Production: Facies, Input, uniform_production
    using CarboKitten.Components.WaterDepth: initial_state
    using CarboKitten.Production: InterpolatedProduction, MultiplyProduction
    <<production-spec>>
end
```
