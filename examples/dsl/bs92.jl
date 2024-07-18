# ~/~ begin <<docs/src/dsl.md#examples/dsl/bs92.jl>>[init]
using CarboKitten.DSL: @spec, @compose

module Units
  using Unitful
  using CarboKitten.Config: Box

  export Amount, Time, Height, Location, Rate, Intensity, Box, @u_str

  const Amount = typeof(1.0u"m")
  const Time = typeof(1.0u"Myr")
  const Height = typeof(1.0u"m")
  const Location = typeof(1.0u"km")
  const Rate = typeof(1.0u"m/Myr")
  const Intensity = typeof(1.0u"W/m^2")
end

# ~/~ begin <<docs/src/dsl.md#dsl-example-time>>[init]
@spec TimeIntegration begin
  using ..Units
  using CarboKitten.Config: TimeProperties

  @kwdef struct Input
    time::TimeProperties
  end

  mutable struct State
    time::Time
  end

  State(input::Input) = State(0.0u"Myr")
end
# ~/~ end
# ~/~ begin <<docs/src/dsl.md#dsl-example-waterdepth>>[init]
@spec WaterDepth begin
  @requires TimeIntegration
  using ..Units

  @kwdef struct Input
    box::Box
    sea_level          # function (t::Time) -> Length
    bedrock_elevation  # function (x::Location, y::Location) -> Length
    subsidence_rate::Rate
  end

  mutable struct State
    sediment_height::Matrix{Height}
  end

  State(input::Input) = State(zeros(Height, input.box.grid_size...))

  function water_depth(input)
    x, y = axes(input.box)
    eta0 = input.bedrock_elevation.(x, y')

    return function(state)
      return input.sea_level(state.time) .- eta0 .+
        (input.subsidence_rate * state.time) .- state.sediment_height
    end
  end
end
# ~/~ end
# ~/~ begin <<docs/src/dsl.md#dsl-example-production>>[init]
@spec UniformProduction begin
  @requires WaterDepth
  using ..Units

  @kwdef struct Facies
    maximum_growth_rate::Rate
    extinction_coefficient::typeof(1.0u"m^-1")
    saturation_intensity::Intensity
  end

  @kwdef struct Input
    insolation::Intensity
    facies::Vector{Facies}
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
end
# ~/~ end

@compose BS92 [UniformProduction]

module BS92
  using CSV
  using DataFrames
  using Interpolations
  using CarboKitten.BoundaryTrait: Shelf

  # In the future, this function will be auto-generated
  function State(input::Input)
    ti_state = TimeIntegration.State(input)
    wd_state = WaterDepth.State(input)
    return State(ti_state.time, wd_state.sediment_height)
  end

  function step(input::Input)
    τ = uniform_production(input)
    function (state::State)
      Δη = τ(state) .* input.time.Δt
      state.sediment_height .+= Δη
      state.time += input.time.Δt
    end
  end

  function sealevel_curve()
       data = DataFrame(CSV.File("data/bs92-sealevel-curve.csv"))
       linear_interpolation(data.time, data.depth)
  end

  const INPUT = Input(
      box = Box(grid_size=(100, 1), phys_scale=600.0u"m"),
      time = TimeProperties(
        Δt = 10u"yr",
        steps = 8000,
        write_interval = 1),
      sea_level = sealevel_curve(),
      bedrock_elevation = (x, y) -> x / 300.0,
      subsidence_rate = 0.0u"m/yr",
      insolation = 400.0u"W/m^2",
      facies = [Facies(
        maximum_growth_rate = 0.005u"m/yr",
        saturation_intensity = 50u"W/m^2",
        extinction_coefficient = 0.05u"m^-1"
      )])

  function run(input::Input)
      step! = step(input)
      state = State(input)
      for i = 1:input.time.steps
          step!(state)
      end
  end
end
# ~/~ end
