# ~/~ begin <<docs/src/transport.md#test/TransportSpec.jl>>[init]
@testset "TransportSpec" begin
  # ~/~ begin <<docs/src/transport.md#transport-spec>>[init]
  using CarboKitten.Config: Box, centered_axes
  using CarboKitten.BoundaryTrait
  using CarboKitten.Transport: deposit, Particle, grid_sample
  using CarboKitten.BoundaryTrait: Periodic

  TestParticle = Particle{Nothing}
  make_vec(x, y) = (x=x, y=y)
  make_test_particle(x, y) = TestParticle(make_vec(x, y), 1.0, 1.0, 1, nothing)
  make_test_particle(x, y, m) = TestParticle(make_vec(x, y), m, 1.0, 1, nothing)

  box = Box{Periodic{2}}(grid_size=(16, 16), phys_scale=1.0 / 16.0 * u"m")
  n_particles = 42
  particles = rand(Float64, 2, n_particles) |> eachcol .|> splat(make_test_particle)

  density = let
    target = zeros(Float64, 1, 16, 16)
    particles .|> deposit(box, target)
    target
  end
  @test sum(density) ≈ n_particles
  # ~/~ end
  # ~/~ begin <<docs/src/transport.md#transport-spec>>[1]
  density = let
    target = zeros(Float64, 1, 16, 16)
    particles = randn(Float64, (2, 1000)) .* 0.1 .+ 0.5 |>
                eachcol .|> splat(make_test_particle)
    particles .|> deposit(box, target)
    target
  end

  let x = ((1:box.grid_size[1]) .- 0.5) .* box.phys_scale |> collect
    mean_x = sum(x .* sum(density; dims=3)[1, :, 1]) ./ sum(density)
    mean_y = sum(x .* sum(density; dims=2)[1, 1, :]) ./ sum(density)
    @test abs(mean_x - 0.5u"m") < 0.01u"m"
    @test abs(mean_y - 0.5u"m") < 0.01u"m"
  end

  let profile(x, y) = exp(-((x - 0.5m)^2 + (y - 0.5m)^2) / 0.04m^2),
      x = ((1:box.grid_size[1]) .- 0.5) .* box.phys_scale |> collect
    density = profile.(centered_axes(box)...)
    particles = [make_test_particle(p.x, p.y, m) for (p, m) in grid_sample(box, density, 1)]

    target = zeros(Float64, 1, 16, 16)
    particles .|> deposit(box, target)

    # test that gaussian blob is centered.
  end
  # ~/~ end
  # ~/~ begin <<docs/src/transport.md#transport-spec>>[2]
  density = let
    target = zeros(Float64, 1, 16, 16)
    Particle((x=0.0, y=0.0), 1.0, 1.0, 1, nothing) |> deposit(box, target)
    target
  end
  @test density[1, 1, 1] ≈ 0.25
  @test density[1, 1, end] ≈ 0.25
  @test density[1, end, 1] ≈ 0.25
  @test density[1, end, end] ≈ 0.25
  # ~/~ end
  # ~/~ begin <<docs/src/transport.md#transport-spec>>[3]
  density = let halfp = (box.phys_scale / 2u"m") |> NoUnits
    target = zeros(Float64, 1, 16, 16)
    Particle((x=halfp, y=halfp),
      1.0, 1.0, 1, nothing) |> deposit(box, target)
    target
  end
  @test density[1, 1, 1] ≈ 1.0
  @test density[1, 1, end] ≈ 0
  @test density[1, end, 1] ≈ 0
  @test density[1, end, end] ≈ 0
  # ~/~ end
  # ~/~ begin <<docs/src/transport.md#transport-spec>>[4]
  height = [0.0 0.0; 1.0 1.0]
  # ~/~ end
end
# ~/~ end
