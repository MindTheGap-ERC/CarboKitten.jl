# Production

```component-dag
CarboKitten.Components.Production
```

The `Production` module specifies the production rate following the model by Bosscher & Schlager 1992 [Bosscher1992](@cite).
The growth rate is given as

$$g(w) = g_m \tanh\left({{I_0 e^{-kw}} \over {I_k}}\right).$$

This can be understood as a smooth transition between the maximum growth rate under saturated conditions, and exponential decay due to light intensity dropping with greater water depth.

``` {.julia #component-production-rate}
function production_rate(insolation, facies, water_depth)
    gₘ = facies.maximum_growth_rate
    I = insolation / facies.saturation_intensity
    x = water_depth * facies.extinction_coefficient
    return water_depth > 0.0u"m" ? gₘ * tanh(I * exp(-x)) : 0.0u"m/Myr"
end
```

From just this equation we can define a uniform production process. This requires that we have a `Facies` that defines the `maximum_growth_rate`, `extinction_coefficient` and `saturation_intensity`.

The `insolation` input may be given as a scalar quantity, say `400u"W/m^2"`, or as a function of time.

``` {.julia file=src/Components/Production.jl}
@compose module Production
@mixin TimeIntegration, WaterDepth, FaciesBase
using ..Common
using ..WaterDepth: water_depth
using ..TimeIntegration: time, write_times
using HDF5
using Logging

export production_rate, uniform_production

@kwdef struct Facies <: AbstractFacies
    maximum_growth_rate::Rate
    extinction_coefficient::typeof(1.0u"m^-1")
    saturation_intensity::Intensity
end

@kwdef struct Input <: AbstractInput
    insolation
end

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

function write_header(fid, input::AbstractInput)
    if input.insolation isa Quantity
        fid["input"]["insolation"] = fill(input.insolation |> in_units_of(u"W/m^2"), input.time.steps + 1)
    elseif input.insolation isa AbstractVector
        fid["input"]["insolation"] = input.insolation |> in_units_of(u"W/m^2")
    else
        t = write_times(input)
        fid["input"]["insolation"] = input.insolation.(t) |> in_units_of(u"W/m^2")
    end

    # fid["input"]["insolation"] = insolation(input) |> in_units_of(u"W/m^2")
    for (i, f) in enumerate(input.facies)
        attr = attributes(fid["input/facies$(i)"])
        attr["maximum_growth_rate"] = f.maximum_growth_rate |> in_units_of(u"m/Myr")
        attr["extinction_coefficient"] = f.extinction_coefficient |> in_units_of(u"m^-1")
        attr["saturation_intensity"] = f.saturation_intensity |> in_units_of(u"W/m^2")
    end
end

<<component-production-rate>>

function uniform_production(input::AbstractInput)
    w = water_depth(input)
    na = [CartesianIndex()]
    insolation_func = insolation(input)
    facies = input.facies

    return function (state::AbstractState)
        return production_rate.(
            insolation_func(state),
            facies[:, na, na],
            w(state)[na, :, :])
    end
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

    function production(input::AbstractInput)
        w = water_depth(input)
        na = [CartesianIndex()]
        output = Array{Amount, 3}(undef, n_facies(input), input.box.grid_size...)

        w = water_depth(input)
        s = insolation(input)
        n_f = n_facies(input)
        facies = input.facies
        Δt = input.time.Δt

        return function(state::AbstractState)
            for f = 1:n_f
                output[f, :, :] = ifelse.(
                    state.ca .== f,
                    production_rate.(s(state), facies, w(state)) .* Δt,
                    0.0u"m")
            end
            return output
        end
    end
end
```

## Tests

``` {.julia file=test/Components/ProductionSpec.jl}
module ProductionSpec
    using Test
    using CarboKitten.Components.Common
    using CarboKitten.Components.Production: Facies, Input, uniform_production
    using CarboKitten.Components.WaterDepth: initial_state

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
end
```
