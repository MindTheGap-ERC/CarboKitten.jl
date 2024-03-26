# ~/~ begin <<docs/src/submarine-transport.md#examples/capt/production-transport.jl>>[init]
module Script
  using CarboKitten
  using CarboKitten.Config: Box, TimeProperties
  using CarboKitten.BoundaryTrait: Shelf, Constant
  using CarboKitten.CaProd: State
  using CarboKitten.CATP: Input, Facies, ProductFrame, submarine_transport, production_propagator
  using CarboKitten.Transport: Box
  using CarboKitten.Burgess2013.CA
  using CarboKitten.Burgess2013: production_rate

  using HDF5
  using Printf
  using ProgressBars
  using Unitful

  using .Iterators: partition, take, map

  const input = Input(
    box = Box{Constant{2,-100}}(grid_size=(100, 50), phys_scale=1.0u"km"),
    time = TimeProperties(
      Δt = 1.0u"kyr",
      steps = 100,
      write_interval = 1
    ),
    sea_level = t -> 4.0u"m" * sin(2π * t / 200.0u"kyr"), 
    subsidence_rate = 50.0u"m/Myr",
    initial_depth = p -> p.x / 2000.0,
    facies = [
        Facies((4, 10), (6, 10), 500.0, 0.800, 300, 1.0, 1.0, 0.05),
        Facies((4, 10), (6, 10), 400.0, 0.1000, 300, 1.0, 1.0, 0.05),
        Facies((4, 10), (6, 10), 100.0, 0.005, 300, 1.0, 1.0, 0.05)
    ],

    insolation = 2000.0,
    Δz = 1.0,
    buffer_depth = 50,
    disintegration_rate = nothing,   # disable disintegration ftm
    wave_shear_stress = nothing,

    g = 9.8,
    transport_subsample = 1,
    transport_max_it = 2000,
    transport_step_size = 0.1,
  )

  function axes(box::Box)
    s = box.phys_scale / u"m" |> NoUnits
    x = collect((0:box.grid_size[1]-1) .* s)
	  y = collect((0:box.grid_size[2]-1) .* s)'
    return x, y
  end

  make_vec2(x, y) = (x=x, y=y)
  phys_grid(box::Box) = make_vec2.(axes(box)...)
  initial_depth(input::Input) = input.box |> phys_grid .|> input.initial_depth

  function updater(input::Input)
    function (s, Δ)
        sr = input.subsidence_rate / u"m/Myr" |> NoUnits
        dt = (input.time.Δt / u"Myr" |> NoUnits)
        dh = (sr .- sum(Δ.production; dims=1)[1,:,:]) .* dt
        s.height .+= dh
        s.time += dt
    end
  end

  function run(input::Input)
    state = State(0.0, initial_depth(input))
    p = production_propagator(input)
    t = submarine_transport(input)
    u = updater(input)

    Channel{ProductFrame}() do ch
      while true
          x = p(state)
          @assert size(x.production)[1] == 3  "production size: $(size(x.production))"
          Δ = t(state, x)
          put!(ch, Δ)
          u(state, Δ)
      end
    end
  end

  function stack_frames(fs::Vector{ProductFrame})  # -> Frame
      ProductFrame(sum(f.production for f in fs))
  end

  function main(input::Input, output::String)
      n_writes = input.time.steps ÷ input.time.write_interval
      x_axis, y_axis = axes(input.box)
      initial_height = initial_depth(input)

      print("Running for $(input.time.steps) time steps and $(n_writes) writes")

      h5open(output, "w") do fid
          gid = create_group(fid, "input")
          gid["x"] = collect(x_axis)
          gid["y"] = collect(y_axis)
          gid["height"] = collect(initial_height)
          dt = input.time.Δt / u"Myr" |> NoUnits
          gid["t"] = collect((0:(n_writes-1)) .* (dt * input.time.write_interval))
          attr = attributes(gid)
          attr["delta_t"] = input.time.Δt / u"Myr" |> NoUnits
          attr["write_interval"] = input.time.write_interval
          attr["time_steps"] = input.time.steps
          attr["subsidence_rate"] = input.subsidence_rate / u"m/Myr" |> NoUnits

          n_facies = length(input.facies)
          ds = create_dataset(fid, "sediment", datatype(Float64),
              dataspace(n_facies, input.box.grid_size..., input.time.steps),
              chunk=(n_facies, input.box.grid_size..., 1))

          results = map(stack_frames, partition(run(input), input.time.write_interval))
          for (step, f) in ProgressBar(enumerate(take(results, n_writes)), total=n_writes)
              ds[:, :, :, step] = f.production
          end
      end
  end
end

Script.main(Script.input, "data/catp.h5")
print("Done\n")
# ~/~ end