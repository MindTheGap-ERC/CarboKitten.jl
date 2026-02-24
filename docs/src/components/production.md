# Production

```component-dag
CarboKitten.Components.Production
```

The `Production` module specifies the production rate following the model by Bosscher & Schlager 1992 [Bosscher1992](@cite).
The growth rate is given as

$$g(w) = g_m \tanh\left({{I_0 e^{-kw}} \over {I_k}}\right).$$

This can be understood as a smooth transition between the maximum growth rate under saturated conditions, and exponential decay due to light intensity dropping with greater water depth.

``` {.julia #component-production-rate}
function benthic_production(insolation, facies, water_depth)
    gₘ = facies.maximum_growth_rate
    I = insolation / facies.saturation_intensity
    x = water_depth * facies.extinction_coefficient
    return x > 0.0 ? gₘ * tanh(I * exp(-x)) : 0.0u"m/Myr"
end

production_rate(i, f, w) = benthic_production(i, f, w)

function capped_production(f, insolation, water_depth, dt)
    clip_positive(x::T) where {T} = max(x, zero(T))
    p = clip_positive(f(insolation, water_depth))
    return min(max(0.0u"m", water_depth), p * dt)
end
```

From just this equation we can define a uniform production process. This requires that we have a `Facies` that defines the `maximum_growth_rate`, `extinction_coefficient` and `saturation_intensity`.

The `insolation` input may be given as a scalar quantity, say `400u"W/m^2"`, or as a function of time.

## Pelagic Production

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
using CarboKitten.Components.Production: benthic_production, pelagic_production
using CarboKitten.ALCAP.Example: INPUT
using CairoMakie

function main()
    water_depth = (0.01:1.0:70.0)u"m"
    fig = Figure()
    facies = [f for f in INPUT.facies if f.active]

    ax_benthic = Axis(fig[1, 1], yreversed=true, title="benthic production", ylabel="depth [m]")
    for f in facies
        p = water_depth .|> (d -> benthic_production(INPUT.insolation, f, d))
        lines!(ax_benthic, p |> in_units_of(u"m/Myr"),
            water_depth |> in_units_of(u"m"), label = string(f.name))
    end

    ax_pelagic = Axis(fig[1, 2], yreversed=true, title="pelagic production")
    for f in facies
        p = water_depth .|> (d -> pelagic_production(INPUT.insolation, f, d))
        lines!(ax_pelagic, p |> in_units_of(u"m^2/Myr"),
            water_depth |> in_units_of(u"m"), label = string(f.name))
    end
    fig[1, 3] = Legend(fig, ax_benthic, "Facies")
    fig
end

end
```

## Production Component

### Production Properties
Because the parameters for benthic and pelagic production have different units, we need different types to store them.

``` {.julia #production-profile}
"""
    production_profile(input::AbstractInput, p)

Given an input production configuration, returns a function of
insolation (in W/m^2) and water depth (in m) to production (in m/Myr).

The default implementation is the identity function; it assumes
that its parameter is a callable object with the correct behaviour.

## Example
If you want to implement your own production profile, define a new type,

    struct MyProduction <: AbstractProduction
        ...
    end

and implement this method for that type,

    import CarboKitten.Components.Production: production_profile

    production_profile(input::AbstractInput, p::MyProduction) =
        function (insolation, water_depth)
            ...
        end
"""
production_profile(input::AbstractInput, p) = p

"""
    is_benthic(obj)

Predicate to determine if a facies or production spec is benthic.
Defaults to `false`.
"""
function is_benthic end

"""
    is_pelagic(obj)

Predicate to determine if a facies or production spec is pelagic.
Defaults to `false`.
"""
function is_pelagic end

struct NoProduction <: AbstractProduction
end

@kwdef struct BenthicProduction <: AbstractProduction
    maximum_growth_rate::Rate = 0.0u"m/Myr"
    extinction_coefficient::typeof(1.0u"m^-1") = 0.0u"m^-1"
    saturation_intensity::Intensity = 1.0u"W/m^2"
end

production_profile(input::AbstractInput, p::BenthicProduction) = (I, w) -> benthic_production(I, p, w)
is_benthic(::BenthicProduction) = true
is_pelagic(::BenthicProduction) = false

@kwdef struct PelagicProduction <: AbstractProduction
    maximum_growth_rate::typeof(1.0u"m/Myr^-1") = 0.0u"Myr^-1"
    extinction_coefficient::typeof(1.0u"m^-1") = 0.0u"m^-1"
    saturation_intensity::Intensity = 1.0u"W/m^2"
    maximum_production_depth::typeof(1.0u"m") = 200.0u"m"
    table_size::Tuple{Int, Int} = (1000, 1000)
end

production_profile(input::AbstractInput, p::PelagicProduction) = (I, w) -> pelagic_production_lookup(I, p, w)
is_benthic(::PelagicProduction) = false
is_pelagic(::PelagicProduction) = true

@kwdef struct Facies <: AbstractFacies
    production = NoProduction()
end

is_benthic(facies::AbstractFacies) = is_benthic(facies.production)
is_pelagic(facies::AbstractFacies) = is_pelagic(facies.production)
```

### Input parameters

In the case of pelagic production, the exact production rate can be quite expensive to compute. Instead, we create look-up tables for each facies. These tables are also easy to store in HDF5. The `maximum_production_depth` parameter controls down to which depth the table is precomputed. At larger depths, we use linear extrapolation (negative values will be clipped to 0). The `production_table_size` parameter determines the resolution of the look-up table.

``` {.julia #production-input}
@kwdef struct Input <: AbstractInput
    insolation
end
```

### Insolation

Insolation can be passed as a constant, function or array (which should have the same length as number of time steps in the run).

``` {.julia #production-insolation}
function insolation(input::AbstractInput)
    insolation = input.insolation
    tprop = input.time
    if insolation isa Quantity
        return s -> insolation
    end
    if insolation isa AbstractVector
        @info "Reading insolation from a table"
        return s -> insolation[s.step+1]
    end
    @info "Reading insolation from a function"
    function (s::AbstractState)
        t = time(tprop, s)
        return insolation(t)
    end
end
```

### Pelagic Lookup tables

We use `linear_interpolation` from `Interpolations` to compute production profiles from look-up tables. If insolation is constant we can use a one-dimensional interpolation. In case of variable insolation, we first find the extrema of the input insolation curve, and then compute a two-dimensional table of values.

``` {.julia #production-lookup}
function pelagic_production_lookup(input::AbstractInput, prod::PelagicProduction)
    insolation_param = input.insolation
    depth_grid = LinRange(0.0, prod.maximum_production_depth, prod.table_size[2])

    if insolation_param isa Quantity
        production_values = [pelagic_production(insolation_param, prod, w) for w in depth_grid]
        interpolated = linear_interpolation(depth_grid, production_values, extrapolation_bc = Line())
        return (I0, w) -> interpolated(w)
    end

    insolation_vec = if insolation_param isa AbstractVector
        insolation_param
    else
        t = time_axis(input)
        insolation_param.(t)
    end
    I_min, I_max = extrema(insolation_vec)

    insolation_grid = LinRange(I_min, I_max, prod.table_size[1])
    production_values = [pelagic_production(I, prod, w) for I in insolation_grid, w in depth_grid]
    return linear_interpolation((insolation_grid, depth_grid), production_values, extrapolation_bc = Line())
end
```

### Boilerplate

``` {.julia file=src/Components/Production.jl}
@compose module Production
@mixin TimeIntegration, WaterDepth, FaciesBase
using ..Common
using ..WaterDepth: water_depth
using ..TimeIntegration: time, write_times
using HDF5
using QuadGK
using Interpolations
using Logging

export production_rate, capped_production, uniform_production, benthic_production, pelagic_production

<<component-production-rate>>
<<pelagic-production>>
<<production-profile>>
<<production-input>>
<<production-insolation>>
<<production-lookup>>

function write_header(input::AbstractInput, output::AbstractOutput)
    if input.insolation isa Quantity
        set_attribute(output, "insolation",
            fill(input.insolation |> in_units_of(u"W/m^2"), input.time.steps))
    elseif input.insolation isa AbstractVector
        set_attribute(output, "insolation",
            input.insolation |> in_units_of(u"W/m^2"))
    else
        t = write_times(input)[1:end-1]
        set_attribute(output, "insolation",
            input.insolation.(t) |> in_units_of(u"W/m^2"))
    end

    for (i, f) in enumerate(input.facies)
        set_attribute(output, "facies$(i)/name", f.name == nothing ? "unnamed" : f.name)
        p = f.production
        if is_pelagic(p)
            set_attribute(output, "facies$(i)/type", "pelagic")
            set_attribute(output, "facies$(i)/maximum_growth_rate", p.maximum_growth_rate |> in_units_of(u"1/Myr"))
            set_attribute(output, "facies$(i)/extinction_coefficient", p.extinction_coefficient |> in_units_of(u"m^-1"))
            set_attribute(output, "facies$(i)/saturation_intensity", p.saturation_intensity |> in_units_of(u"W/m^2"))
        elseif is_benthic(p)
            set_attribute(output, "facies$(i)/type", "benthic")
            set_attribute(output, "facies$(i)/maximum_growth_rate", p.maximum_growth_rate |> in_units_of(u"m/Myr"))
            set_attribute(output, "facies$(i)/extinction_coefficient", p.extinction_coefficient |> in_units_of(u"m^-1"))
            set_attribute(output, "facies$(i)/saturation_intensity", p.saturation_intensity |> in_units_of(u"W/m^2"))
        end
    end
end

function uniform_production(input::AbstractInput)
    w = water_depth(input)
    na = [CartesianIndex()]
    insolation_func = insolation(input)
    facies = input.facies
    dt = input.time.Δt
    production_rates = [production_profile(input, f.production) for f in facies]

    p(state::AbstractState, wd::AbstractMatrix) =
        capped_production.(production_rates[:, na, na], insolation_func(state), wd[na, :, :], dt)

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

``` {.julia file=src/Components/CAProduction.jl}
@compose module CAProduction
    @mixin TimeIntegration, CellularAutomaton, Production
    using ..Common
    using ..Production: production_rate, insolation
    using ..WaterDepth: water_depth
    using Logging

    function production(input::AbstractInput)
        w = water_depth(input)
        na = [CartesianIndex()]
        output_ = Array{Amount, 3}(undef, n_facies(input), input.box.grid_size...)

        w = water_depth(input)
        s = insolation(input)
        n_f = n_facies(input)
        facies = input.facies
        active_facies = [i for (i, f) in pairs(facies) if f.active]
        global_facies = [i for (i, f) in pairs(facies) if !f.active]
        dt = input.time.Δt
        # Having this a Tuple should make things type stable?
        production_specs = ((production_profile(input, f.production) for f in facies)...,)

        function p(state::AbstractState, wd::AbstractMatrix)::Array{Amount,3}
            output::Array{Amount, 3} = output_
            insolation::typeof(1.0u"W/m^2") = s(state)
            for i in eachindex(IndexCartesian(), wd)
                for f in eachindex(facies)
                    if facies[f].active
                        output[f, i[1], i[2]] = f != state.ca[i] ? 0.0u"m" :
                            capped_production(production_specs[f], insolation, wd[i], dt)
                    else
                        output[f, i[1], i[2]] =
                            capped_production(production_specs[f], insolation, wd[i], dt)
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

### If production is higher in shallower water 

```{.julia #production-spec}
@testset "Components/Production" begin
    let facies = Facies(
            maximum_growth_rate = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        input = Input(
            box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
            time = TimeProperties(Δt=1.0u"kyr", steps=10),
            sea_level = t -> 0.0u"m",
            initial_topography = (x, y) -> -10u"m",
            subsidence_rate = 0.0u"m/Myr",
            facies = [facies],
            insolation = 400.0u"W/m^2")

        state = initial_state(input)
        prod = uniform_production(input)(state)
        @test all(prod[1:end-1,:] .>= prod[2:end,:])
    end
end
```

### Variable insolation

If insolation increases linearly with time, production at t = 10 should be higher than at t = 1.

```{.julia #production-spec}
@testset "Components/Production/variable_insolation" begin
    let facies = Facies(
            maximum_growth_rate = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        input = Input(
            box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
            time = TimeProperties(Δt=1.0u"kyr", steps=10),
            sea_level = t -> 0.0u"m",
            initial_topography = (x, y) -> -10u"m",
            subsidence_rate = 0.0u"m/Myr",
            facies = [facies],
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
    using CarboKitten.Components.Common
    using CarboKitten.Components.Production: Facies, Input, uniform_production
    using CarboKitten.Components.WaterDepth: initial_state

    <<production-spec>>
end
```
