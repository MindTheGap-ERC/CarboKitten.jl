# ~/~ begin <<docs/src/output/abstract.md#test/Output/WaterDepthSpec.jl>>[init]
using CarboKitten
using CarboKitten.Models: WithoutCA as M
using CarboKitten.Output: bedrock
using CarboKitten.Export: read_volume

function test_model_input(;
        box = CarboKitten.Box{Coast}(grid_size=(100, 1), phys_scale=150.0u"m"),
        subsidence_rate = 38.0u"m/Myr")
    facies = [
        M.Facies(
            production=CarboKitten.Production.EXAMPLE[:euphotic] * ProductionBoost(factor=0.25),
            transport_coefficient=50.0u"m/Myr"),
        M.Facies(
            production=CarboKitten.Production.EXAMPLE[:oligophotic] * ProductionBoost(factor=0.25),
            transport_coefficient=25.0u"m/Myr"),
        M.Facies(
            production=CarboKitten.Production.EXAMPLE[:aphotic] * ProductionBoost(factor=0.25),
            transport_coefficient=12.5u"m/Myr"),
    ]

    M.Input(
        box = box,
        time = TimeProperties(Δt=0.0002u"Myr", steps=5000),
        facies = facies,
        initial_topography = (x, y) -> -x / 300.0,

        sea_level = t -> 5.5u"m" * sin(2π * t / 0.18u"Myr"),
        subsidence_rate = subsidence_rate,
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

@testset "Output.Abstract.water_depth" begin
    na = [CartesianIndex()]
    output = test_model()
    section = output.data_volumes[:full][:, 1]
    wd = water_depth(output.header, section)
    sc = stratigraphic_column(section)
    st_preserved = dropdims(cumsum(sum(sc; dims=1); dims=3), dims=1)
    st = dropdims(
        cumsum(sum(section.deposition .- section.disintegration; dims=1), dims=3),
        dims=1)

    @test st_preserved[:, end] ≈ st[:, end]

    initial_bathymetry = output.header.initial_topography[:, na]
    subsidence = output.header.subsidence_rate .* output.header.axes.t[na, :]
    reconstructed_bathymetry = initial_bathymetry .+ st .- subsidence

    @test reconstructed_bathymetry ≈ section.bathymetry

    sea_level = output.header.sea_level[na, :]
    reconstructed_wd = sea_level .- reconstructed_bathymetry

    @test reconstructed_wd ≈ wd
end
# ~/~ end
# ~/~ begin <<docs/src/output/abstract.md#test/Output/WaterDepthSpec.jl>>[1]
@testset "Output.Abstract.variable_subsidence" begin
    σ(x, y, t) = (20.0u"m/Myr" + x * 5.0u"1/Myr" + y * 10.0u"1/Myr") *
        exp(-t / 0.5u"Myr")
    input = test_model_input(
        box = CarboKitten.Box{Coast}(grid_size=(50, 3), phys_scale=150.0u"m"),
        subsidence_rate = SpaceTimeVariableSubsidenceRate(σ))

    x, y = box_axes(input.box)
    t = time_axis(input.time)
    rates = σ.(x, y', reshape(t[1:end-1], 1, 1, :))
    subsidence = cat(zeros(typeof(1.0u"m"), size(x)..., size(y)..., 1),
                     cumsum(rates .* input.time.Δt; dims=3); dims=3)

    function check(header, data)
        @test header.subsidence_rate === nothing
        @test header.total_subsidence ≈ subsidence[:, :, end]
        @test bedrock(header) ≈ header.initial_topography .- subsidence[:, :, end]

        st = dropdims(
            cumsum(sum(data.deposition .- data.disintegration; dims=1), dims=4),
            dims=1)
        reconstructed_bathymetry = header.initial_topography .+ st .- subsidence
        @test reconstructed_bathymetry ≈ data.bathymetry
    end

    output = MemoryOutput(input)
    run_model(Model{M}, input, output)
    check(output.header, output.data_volumes[:full])

    mktempdir() do path
        filename = joinpath(path, "variable_subsidence.h5")
        run_model(Model{M}, input, filename)
        header, data = read_volume(filename, :full)
        check(header, data)
    end
end
# ~/~ end
