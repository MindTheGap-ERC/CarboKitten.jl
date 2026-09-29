# ~/~ begin <<docs/src/visualization/profiles.md#test/Visualization/UnconformitiesSpec.jl>>[init]
using CarboKitten
using CarboKitten.Models: WithoutCA as M

function test_model_input()
    facies = [
        M.Facies(
            production=CarboKitten.Production.EXAMPLE[:euphotic],
            transport_coefficient=50.0u"m/Myr"),
        M.Facies(
            production=CarboKitten.Production.EXAMPLE[:oligophotic],
            transport_coefficient=25.0u"m/Myr"),
        M.Facies(
            production=CarboKitten.Production.EXAMPLE[:aphotic],
            transport_coefficient=12.5u"m/Myr"),
    ]

    M.Input(
        box = CarboKitten.Box{Coast}(grid_size=(100, 1), phys_scale=150.0u"m"),
        time = TimeProperties(Δt=0.0002u"Myr", steps=5000),
        facies = facies,
        initial_topography = (x, y) -> -x / 300.0,

        sea_level = t -> 5.5u"m" * sin(2π * t / 0.18u"Myr"),
        subsidence_rate = 38.0u"m/Myr",
        disintegration_rate = 80.0u"m/Myr",
        lithification_time = 100.0u"yr",
        insolation = 400.0u"W/m^2",
    )
end

function test_model()
    input = test_model_input()
    output = MemoryOutput(input)
    run_model(Model{M}, input, output)
end

@testset "CarboKitten.Visualization" begin
end
# ~/~ end
