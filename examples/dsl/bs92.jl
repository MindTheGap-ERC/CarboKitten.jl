# ~/~ begin <<docs/src/dsl.md#examples/dsl/bs92.jl>>[init]
using ModuleMixins: @compose

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

  export AbstractInput, AbstractFacies, AbstractState

  abstract type AbstractInput end
  abstract type AbstractFacies end
  abstract type AbstractState end
end

# ~/~ begin <<docs/src/dsl.md#dsl-example-time>>[init]
@compose module TimeIntegration
  using ..Units
  using CarboKitten.Config: TimeProperties

  @kwdef struct Input <: AbstractInput
    time::TimeProperties
  end

  mutable struct State <: AbstractState
    time::Time
  end

  State(input::AbstractInput) = State(0.0u"Myr")
end
# ~/~ end
# ~/~ begin <<docs/src/dsl.md#dsl-example-waterdepth>>[init]
@compose module WaterDepth
  @mixin TimeIntegration
  using ..Units
  using CarboKitten.Config: axes

  @kwdef struct Input <: AbstractInput
    box::Box
    sea_level          # function (t::Time) -> Length
    bedrock_elevation  # function (x::Location, y::Location) -> Length
    subsidence_rate::Rate
  end

  mutable struct State <: AbstractState
    sediment_height::Matrix{Height}
  end

  State(input::AbstractInput) = State(zeros(Height, input.box.grid_size...))

  function water_depth(input::AbstractInput)
    x, y = axes(input.box)
    eta0 = input.bedrock_elevation.(x, y')

    return function(state::AbstractState)
      return input.sea_level(state.time) .- eta0 .+
        (input.subsidence_rate * state.time) .- state.sediment_height
    end
  end
end
# ~/~ end
# ~/~ begin <<docs/src/dsl.md#dsl-example-production>>[init]
@compose module UniformProduction
  @mixin WaterDepth
  using ..Units
  using ..WaterDepth: water_depth

  @kwdef struct Facies <: AbstractFacies
    maximum_growth_rate::Rate
    extinction_coefficient::typeof(1.0u"m^-1")
    saturation_intensity::Intensity
  end

  @kwdef struct Input <: AbstractInput
    insolation::Intensity
    facies::Vector{Facies}
  end

  function production_rate(insolation, facies, water_depth)
      gₘ = facies.maximum_growth_rate
      I = insolation / facies.saturation_intensity
      x = water_depth * facies.extinction_coefficient
      return water_depth > 0.0u"m" ? gₘ * tanh(I * exp(-x)) : 0.0u"m/Myr"
  end

  function uniform_production(input::AbstractInput)
    w = water_depth(input)
    na = [CartesianIndex()]

    return function(state::AbstractState)
      return production_rate.(
        input.insolation,
        input.facies[:,na,na],
        w(state)[na,:,:])
    end
  end
end
# ~/~ end

@compose module BS92
  @mixin UniformProduction

  using CSV
  using DataFrames
  using Interpolations
  using CarboKitten.BoundaryTrait: Shelf
  using ..UniformProduction: uniform_production
  using ..TimeIntegration
  using ..WaterDepth

  # In the future, this function will be auto-generated
  function State(input::Input)
    ti_state = TimeIntegration.State(input)
    wd_state = WaterDepth.State(input)
    return State(ti_state.time, wd_state.sediment_height)
  end

  function step(input::Input)
    τ = uniform_production(input)
    function (state::State)
      Δη = sum(τ(state); dims=1)[1,:,:] .* input.time.Δt
      state.sediment_height .+= Δη
      state.time += input.time.Δt
    end
  end

  function sealevel_curve()
       data = DataFrame(CSV.File("data/bs92-sealevel-curve.csv"))
       linear_interpolation(data.time, data.depth)
  end

  const INPUT = Input(
      box = Box{Shelf}(grid_size=(100, 1), phys_scale=600.0u"m"),
      time = TimeProperties(
        Δt = 10.0u"yr",
        steps = 8000,
        write_interval = 100),
      sea_level = let sc = sealevel_curve()
        t -> -sc(t / u"yr") * u"m"
      end,
      bedrock_elevation = (x, y) -> - x / 300.0,
      subsidence_rate = 0.0u"m/yr",
      insolation = 400.0u"W/m^2",
      facies = [Facies(
        maximum_growth_rate = 0.005u"m/yr",
        saturation_intensity = 50.0u"W/m^2",
        extinction_coefficient = 0.05u"m^-1"
      )])

  function run(input::Input)
      step! = step(input)
      getwd = WaterDepth.water_depth(input)
      state = State(input)

      n_writes = input.time.steps ÷ input.time.write_interval
      result = Array{Amount, 2}(undef, input.box.grid_size[1], n_writes)
      wd = Array{Amount, 2}(undef, input.box.grid_size[1], n_writes)
      for i = 1:n_writes
        wd[:,i] = getwd(state)
        for _ = 1:input.time.write_interval
              step!(state)
          end
          result[:,i] = state.sediment_height[:,1]
      end
      return result, wd
  end
end

using GLMakie
using CarboKitten.Utility: in_units_of
using CarboKitten.Config: axes as box_axes
using Unitful

function main()
  result, _ = BS92.run(BS92.INPUT)
  fig = Figure()
  ax = Axis(fig[1,1], xlabel="x (km)", ylabel="z (m)")
  x, y = box_axes(BS92.INPUT.box)
  η0 = BS92.INPUT.bedrock_elevation.(x, y')

  for l in eachcol(result)
    η = η0 .+ l
    lines!(ax, x |> in_units_of(u"km"), vec(η) |> in_units_of(u"m"), color=:steelblue4)
  end

  fig
end

main()
# ~/~ end
