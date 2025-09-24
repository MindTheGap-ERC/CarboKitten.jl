using CarboKitten
using Test
using Unitful

@testset "Unit test sanity" begin
    @test 1 + 1 == 2
end

include("Unitful.jl")

TEST_PATH = mktempdir()

@testset "CarboKitten" begin
    CarboKitten.init()
    @info "Running models for property testing"

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
          maximum_growth_rate = 0.005u"m/yr",
          saturation_intensity = 50.0u"W/m^2",
          extinction_coefficient = 0.05u"m^-1"
        )])

    run_model(Model{BS92}, BS92_TEST_INPUT, joinpath(TEST_PATH, "bs92_spm.h5"))

    include("UtilitySpec.jl")
    include("StencilSpec.jl")
    include("DenudationSpec.jl")
    include("SedimentStackSpec.jl")

    include("Components/CellularAutomatonSpec.jl")
    include("Components/FaciesBaseSpec.jl")
    include("Components/TimeIntegrationSpec.jl")
    include("Components/BoxesSpec.jl")
    include("Components/ProductionSpec.jl")
    include("Components/InitialSedimentSpec.jl")
    include("Components/H5WriterSpec.jl")

    include("Transport/AdvectionSpec.jl")
    include("Transport/IntertidalZoneSpec.jl")

    include("ExportSpec.jl")
    include("OutputDataSpec.jl")
end
