# ~/~ begin <<docs/src/dsl.md#examples/dsl/bs92.jl>>[init]
using CarboKitten.DSL: @spec, @compose

module Units
  using Unitful
  using CarboKitten.Config: Box

  export Amount, Time, Height, Location, Rate, Intensity, Box

  const Amount = typeof(1.0u"m")
  const Time = typeof(1.0u"Myr")
  const Height = typeof(1.0u"m")
  const Location = typeof(1.0u"km")
  const Rate = typeof(1.0u"m/Myr")
  const Intensity = typeof(1.0u"W/m^2")
end

# ~/~ begin <<docs/src/dsl.md#dsl-example-time>>[init]
@spec TimeIntegration begin
  using ..Types
  using CarboKitten.Config: TimeProperties

  struct Input
    time::TimeProperties
  end

  struct State
    time::Time
  end
end
# ~/~ end
# ~/~ begin <<docs/src/dsl.md#dsl-example-waterdepth>>[init]
@spec WaterDepth begin
  @requires TimeIntegration
  using ..Types

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

function water_depth(input)
  x, y = axes(input.box)
  eta0 = input.bedrock_elevation.(x, y')

  return function(state::State)
    return input.sea_level(state.time) .- eta0 .+
      (input.subsidence_rate * state.time) .- state.sediment_height
  end
end
# ~/~ end
# ~/~ begin <<docs/src/dsl.md#dsl-example-production>>[init]
@spec UniformProduction begin
  @requires WaterDepth
  using ..Types

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

function uniform_production(input)
  w = water_depth(input)
  na = [CartesianIndex()]

  return function(state)
    return production_rate.(
      input.insolation,
      input.facies[:,na,na],
      w(state)[na,:,:])
  end
end
# ~/~ end

@compose BS92 [UniformProduction]

module BS92

end
# ~/~ end