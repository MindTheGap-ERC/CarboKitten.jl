# ~/~ begin <<docs/src/submarine-transport.md#examples/capt/production-transport.jl>>[init]
module Script
  using CarboKitten
  using CarboKitten.BoundaryTrait: Shelf
  using CarboKitten.CaProd: State
  using CarboKitten.CATP: Input, Facies
  using CarboKitten.Transport: Box

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
        Facies((4, 10), (6, 10), 500.0, 0.8, 300, 1.0, 1.0, 1.0),
        Facies((4, 10), (6, 10), 400.0, 0.1, 300, 1.0, 1.0, 1.0),
        Facies((4, 10), (6, 10), 100.0, 0.005, 300, 1.0, 1.0, 1.0)
    ],

    insolation = 2000.0,
    Δz = 1.0,
    buffer_depth = 50,
    disintegration_rate = Nothing,   # disable disintegration ftm
    wave_shear_stress = Nothing,

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

  function run(input::Input)
    box = make_box(input)
    state = State()
    Channel{Frame}() do ch
      while true
          Δ = p(state)
          put!(ch, Δ)
          u(state, Δ)
      end
    end
  end
end

# Script.run()
# ~/~ end