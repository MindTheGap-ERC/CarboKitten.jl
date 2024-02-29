# ~/~ begin <<docs/src/submarine-transport.md#examples/capt/production-transport.jl>>[init]
module Script
  using CarboKitten
  using CarboKitten.BoundaryTrait: Shelf
  using CarboKitten.CaProd: State
  using CarboKitten.CATP: Input, Facies, ProductFrame, submarine_transport, production_propagator
  using CarboKitten.Transport: Box
  using CarboKitten.Burgess2013.CA
  using CarboKitten.Burgess2013: production_rate

  using HDF5
  using Printf

  using .Iterators: partition, take

  const input = Input(
    sea_level = t -> 4 * sin(2π * t / 0.2), 
    subsidence_rate = 50.0,
    initial_depth = p -> p.x / 2,
    grid_size = (100, 50),
    boundary = Shelf,
    phys_scale = 1.0,
    Δt = 0.0001,
    time_steps = 100,
    write_interval = 1,
    facies = [
        Facies((4, 10), (6, 10), 500.0, 0.8, 300, 1.0, 1.0, 10.0),
        Facies((4, 10), (6, 10), 400.0, 0.1, 300, 1.0, 1.0, 10.0),
        Facies((4, 10), (6, 10), 100.0, 0.005, 300, 1.0, 1.0, 10.0)
    ],

    insolation = 2000.0,
    Δz = 1.0,
    buffer_depth = 50,
    disintegration_rate = nothing,   # disable disintegration ftm
    wave_shear_stress = nothing,

    g = 9.8,
    transport_subsample = 1 
  )

  function make_box(input)
    phys_size = (x = input.grid_size[1] * input.phys_scale,
                 y = input.grid_size[2] * input.phys_scale)
    Box(input.grid_size, phys_size, input.phys_scale)
  end

  function axes(box::Box)
    x = collect((0:box.grid_size[1]-1) .* box.phys_scale)
	  y = collect((0:box.grid_size[2]-1) .* box.phys_scale)'
    return x, y
  end

  make_vec2(x, y) = (x=x, y=y)
  phys_grid(box::Box) = make_vec2.(axes(box)...)
  initial_depth(input::Input) = input |> make_box |> phys_grid .|> input.initial_depth

  function updater(input::Input)
    function (s, Δ)
        s.height .+= (input.subsidence_rate .- sum(Δ.production; dims=1)[1,:,:]) * input.Δt
        s.time += input.Δt
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
      n_writes = input.time_steps ÷ input.write_interval
      x_axis, y_axis = axes(make_box(input))
      initial_height = initial_depth(input)

      h5open(output, "w") do fid
          gid = create_group(fid, "input")
          gid["x"] = collect(x_axis)
          gid["y"] = collect(y_axis)
          gid["height"] = collect(initial_height)
          gid["t"] = collect((0:(n_writes-1)) .* (input.Δt * input.write_interval))
          attr = attributes(gid)
          attr["delta_t"] = input.Δt
          attr["write_interval"] = input.write_interval
          attr["time_steps"] = input.time_steps
          attr["subsidence_rate"] = input.subsidence_rate

          n_facies = length(input.facies)
          ds = create_dataset(fid, "sediment", datatype(Float64),
              dataspace(n_facies, input.grid_size..., input.time_steps),
              chunk=(n_facies, input.grid_size..., 1))

          results = map(stack_frames, partition(run(input), input.write_interval))
          for (step, f) in enumerate(take(results, n_writes))
              ds[:, :, :, step] = f.production
          end
      end
  end
end

Script.main(Script.input, "data/catp.h5")
print("Done\n")
# ~/~ end