# The CarboKitten DSL

At some point in the evolution of CarboKitten, there will be many different implementations for production, disintegration, transport etc. Each with slightly different needs on the `Input` and `State` structures. Wouldn't it be cool if we could compose a model directly out of these components?

``` {.julia file=src/DSL.jl}
module DSL

using Unitful

<<dsl>>

end
```

## Some types

``` {.julia #dsl}
const Amount = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")
const Height = typeof(1.0u"m")
const Location = typeof(1.0u"km")
const Rate = typeof(1.0u"m/Myr")
const Intensity = typeof(1.0u"W/m^2")

abstract type Input end
abstract type State end
abstract type Facies end
```

## Water depth

``` {.julia #dsl}
@spec WaterDepth begin
  struct Input
    box::Box
    sea_level          # function (t::Time) -> Length
    bedrock_elevation  # function (x::Location, y::Location) -> Length
    subsidence_rate::Rate
  end

  struct State
    time::Time
    sediment_height::Matrix{Height}
  end
end

function water_depth(input::Input)
  x, y = axes(input.box)
  eta0 = input.bedrock_elevation.(x, y')

  return function(state::State)
    return input.sea_level(state.time) .- eta0 .+
      (input.subsidence_rate * state.time) .- state.sediment_height
  end
end
```

## Uniform Production

``` {.julia #dsl}
@spec UniformProduction begin
  @requires WaterDepth

  struct Facies
    maximum_growth_rate::Rate
    extinction_coefficient::typeof(1.0u"m^-1")
    saturation_intensity::Intensity
  end

  struct Input
    insolation::Intensity
    facies::Vector{Facies}
  end
end

function production_rate(insolation, facies, water_depth)
    gₘ = facies.maximum_growth_rate
    I = insolation / facies.saturation_intensity
    x = water_depth * facies.extinction_coefficient
    return water_depth > 0.0u"m" ? gₘ * tanh(I * exp(-x)) : 0.0u"m/Myr"
end

function uniform_production(input::Input)
  w = water_depth(input)
  na = [CartesianIndex()]

  return function(state::State)
    return production_rate.(
      input.insolation,
      input.facies[:,na,na],
      w(state)[na,:,:])
  end
end
```
