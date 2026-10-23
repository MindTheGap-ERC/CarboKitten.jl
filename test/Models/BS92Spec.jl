# ~/~ begin <<docs/src/test-cases.md#test/Models/BS92Spec.jl>>[init]
using CarboKitten

@testset "Cases/BS92" begin
    TEST_PATH = mktempdir()
    BS92_TEST_INPUT = BS92.Input(
        tag = "single pixel model",
        box = Box{Periodic{2}}(grid_size=(1, 1), phys_scale=600.0u"m"),
        time = TimeProperties(
          Δt = 10.0u"yr",
          steps = 8000),
        output = Dict(:full => OutputSpec((1, 1), 80)),
        sea_level = t -> 10.0u"m" * sin(2π * t / 20u"kyr"),
        initial_topography = (_, _) -> - 50.0u"m",
        subsidence_rate = 0.001u"m/yr",
        insolation = 400.0u"W/m^2",
        facies = [BS92.Facies(
            production = BenthicProduction(
                maximum_growth_rate = 0.005u"m/yr",
                saturation_intensity = 50.0u"W/m^2",
                extinction_coefficient = 0.05u"m^-1"))
        ])
    run_model(Model{BS92}, BS92_TEST_INPUT, joinpath(TEST_PATH, "bs92_spm.h5"))

    @testset "Water depth signs" begin
        using CarboKitten.Export: extract_wd, extract_sac, read_column

        test_path::String = TEST_PATH
        header, data = read_column(joinpath(test_path, "bs92_spm.h5"), :full)
        wd = extract_wd(header, data, 1)
        sac = extract_sac(header, data, 1)
        submerged = wd.wd_1 .> -1.0u"m"
        growing = (sac.sac_1[2:end] .- sac.sac_1[1:end-1]) .> 0.5u"m"
        @test all(growing .&& (submerged[1:end-1] .|| submerged[2:end]) .|| .!growing)
    end
end
# ~/~ end
