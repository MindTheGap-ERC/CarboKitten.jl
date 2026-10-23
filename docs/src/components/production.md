# Production

```component-dag
CarboKitten.Components.Production
```

## Production Component

### Production Properties

### Input parameters

``` {.julia #production-input}
@kwdef struct Input <: AbstractInput
    insolation = 400.0u"W/m^2"
end

@kwdef struct Facies <: AbstractFacies
    production = NoProduction()
end

is_benthic(facies::AbstractFacies) = is_benthic(facies.production)
is_pelagic(facies::AbstractFacies) = is_pelagic(facies.production)
is_interpolated(facies::AbstractFacies) = is_interpolated(facies.production)
```

The `production_modifiers` field has been removed from `Input`. Time-dependent
behaviour is now expressed by composing production specs directly in the `Facies`
definition using `MultiplyProduction`.

## Production cap

Because we can only produce as much as keeps the factory submerged, we have to cap the total production in a single time step to the current water depth. This assumes we have production as a function of time and water depth.

``` {.julia #capped-production}
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

``` {.julia file=src/Components/Production.jl}
@compose module Production
@mixin TimeIntegration, WaterDepth, FaciesBase
using ..Common
using ..WaterDepth: water_depth
using ..TimeIntegration: time, write_times
using ...Production: NoProduction, production_profile, insolation_curve
import ...Production: is_benthic, is_pelagic, is_interpolated

using HDF5
using QuadGK
using Interpolations
using Logging

export uniform_production

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

<<capped-production>>

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
    using ...Production: production_profile
    using Logging

    function production(input::AbstractInput)
        w = water_depth(input)
        na = [CartesianIndex()]
        output_ = Array{Amount, 3}(undef, n_facies(input), input.box.grid_size...)

        facies = input.facies
        dt = input.time.Δt
        production_specs = ((Production.production_profile(input, f.production) for f in facies)...,)
        get_time = time(input)

        function p(state::AbstractState, wd::AbstractMatrix)::Array{Amount,3}
            output::Array{Amount, 3} = output_
            t = get_time(state)
            for i in eachindex(IndexCartesian(), wd)
                for f in eachindex(facies)
                    if facies[f].active
                        output[f, i[1], i[2]] = f != state.ca[i] ? 0.0u"m" :
                            Production.capped_production(production_specs[f], t, wd[i], dt)
                    else
                        output[f, i[1], i[2]] =
                            Production.capped_production(production_specs[f], t, wd[i], dt)
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
