# ~/~ begin <<docs/src/transport.md#test/TransportSpec.jl>>[init]
@testset "TransportSpec" begin
    # ~/~ begin <<docs/src/transport.md#transport-spec>>[init]
    using CarboKitten.Transport: deposit, Box, Particle
    using CarboKitten.BoundaryTrait: Periodic

    TestParticle = Particle{Nothing}

    box = Box((16, 16), (x=1.0, y=1.0), 1.0/16.0)
    n_particles = 42
    particles = rand(Float64, 2, n_particles) |> eachcol .|> 
    		(p -> TestParticle((x=p[1], y=p[2]), 1.0, 1.0, 1, nothing));

    density = let
    	target = zeros(Float64, 1, 16, 16)
    	particles .|> deposit(Periodic{2}, box, target)
    	target
    end
    @test sum(density) ≈ n_particles
    # ~/~ end
    # ~/~ begin <<docs/src/transport.md#transport-spec>>[1]
    density = let
        target = zeros(Float64, 1, 16, 16)
        Particle((x=0.0, y=0.0), 1.0, 1.0, 1, nothing) |> deposit(Periodic{2}, box, target)
        target
    end
    @test density[1,1,1] ≈ 0.25
    @test density[1,1,end] ≈ 0.25
    @test density[1,end,1] ≈ 0.25
    @test density[1,end,end] ≈ 0.25
    # ~/~ end
    # ~/~ begin <<docs/src/transport.md#transport-spec>>[2]
    density = let
        target = zeros(Float64, 1, 16, 16)
        Particle((x=box.phys_scale/2, y=box.phys_scale/2),
            1.0, 1.0, 1, nothing) |> deposit(Periodic{2}, box, target)
        target
    end
    @test density[1,1,1] ≈ 1.0
    @test density[1,1,end] ≈ 0
    @test density[1,end,1] ≈ 0
    @test density[1,end,end] ≈ 0
    # ~/~ end
    # ~/~ begin <<docs/src/transport.md#transport-spec>>[3]
    height = [0.0 0.0; 1.0 1.0]
    # ~/~ end
end
# ~/~ end