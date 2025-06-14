# ~/~ begin <<docs/src/data-export.md#test/ExportSpec.jl>>[init]
using CarboKitten
using CarboKitten.Export: Axes, Header, DataVolume, data_export, CSVExportTrait,
    age_depth_model, extract_sac, extract_sc, CSV, read_data, extract_sac, extract_wd
using CSV: read as read_csv
using TOML
using DataFrames
using Unitful

const Amount = typeof(1.0u"m")

# ~/~ begin <<docs/src/data-export.md#export-test-case>>[init]
const AXES1 = Axes(
    x=[0.0, 1.0, 2.0] * u"m",
    y=[1.0] * u"m",
    t=(0.0:0.1:1.0) * u"Myr")

const HEADER1 = Header(
    tag="test",
    axes=AXES1,
    Δt=0.1u"Myr",
    time_steps=10,
    initial_topography=zeros(typeof(1.0u"m"), 3, 3),
    sea_level=zeros(typeof(1.0u"m"), 11),
    subsidence_rate=10u"m/Myr")

const PRODUCTION1 = reshape(
    hcat(ones(Amount, 10),
        ones(Amount, 10),
        cumsum(ones(Amount, 10)) / 5.5)',
    1, 3, 1, 10)

const DISINTEGRATION1 = reshape(
    hcat(zeros(Amount, 10),
        1:10 .|> (x -> x < 4 || x > 6 ? 0.0u"m" : 2.0u"m"),
        zeros(Amount, 10))',
    1, 3, 1, 10)

const ELEVATION1 = cat(
    [0.0, 0.0, 0.0]u"m",
    cumsum(PRODUCTION1 .- DISINTEGRATION1; dims=4)[1, :, :, :];
    dims=3)

const DATA1 = DataVolume(
    slice=(:,:),
    write_interval=1,
    disintegration=DISINTEGRATION1,
    production=PRODUCTION1,
    deposition=PRODUCTION1 .- DISINTEGRATION1,
    sediment_thickness=ELEVATION1)

const GRID_LOCATIONS1 = [(1, 1), (2, 1), (3, 1)]

const COLUMNS1 = [DATA1[loc...] for loc in GRID_LOCATIONS1]
# ~/~ end

@testset "Data Export" begin
    # ~/~ begin <<docs/src/data-export.md#export-test>>[init]
    @testset "Hither and Dither" begin
        df = data_export(CSVExportTrait{:sediment_accumulation_curve},
                         HEADER1, COLUMNS1)
        @test df.sac_1 ≈ ELEVATION1[1, 1, :]
        @test df.sac_2 ≈ ELEVATION1[2, 1, :]
        @test df.sac_3 ≈ ELEVATION1[3, 1, :]
    end
    # ~/~ end
    # ~/~ begin <<docs/src/data-export.md#export-test>>[1]
    @testset "ADM Monotonicity" begin
        sac = data_export(CSVExportTrait{:sediment_accumulation_curve}, HEADER1, COLUMNS1)
        adm = data_export(CSVExportTrait{:age_depth_model}, HEADER1, COLUMNS1)

        @test sac.sac_1 == adm.adm_1
        @test sac.sac_3 == adm.adm_3
        @test sac.sac_2 != adm.adm_2

        @test all(adm.adm_2[2:end] .- adm.adm_2[1:end-1] .>= 0.0u"m")
    end
    # ~/~ end
    # ~/~ begin <<docs/src/data-export.md#export-test>>[2]
    # @testset "SC sum equals ADM" begin
    #     sac = extract_sac(HEADER1, DATA1, COLUMNS1)
    #     adm = sac |> age_depth_model
    #     sc = extract_sc(HEADER1, DATA1, COLUMNS1)
    #     @test [0.0u"m"; cumsum(sc.sc1_f1)] ≈ adm.adm1
    #     @test [0.0u"m"; cumsum(sc.sc2_f1)] ≈ adm.adm2
    #     @test [0.0u"m"; cumsum(sc.sc3_f1)] ≈ adm.adm3
    # end
    # ~/~ end
    # @testset "Write to folder" begin
    #     mktempdir() do path
    #         spec = CSV(
    #             :sediment_accumulation_curve => joinpath(path, "sac.csv"),
    #             :age_depth_model => joinpath(path, "adm.csv"),
    #             :stratigraphic_column => joinpath(path, "sc.csv"),
    #             :water_depth => joinpath(path, "wd.csv"),
    #             :metadata => joinpath(path, "metadata.toml"))
    #         data_export(spec, HEADER1, COLUMNS1)
    #         for f in values(spec.output_files)
    #             @test isfile(f)
    #         end

    #         metadata = TOML.parsefile(spec.output_files[:metadata])
    #         @test IdDict(Symbol(k) => v for (k, v) in metadata["files"]) == spec.output_files
    #         @test length(metadata["locations"]) == 3
    #         adm = read_csv(spec.output_files[:age_depth_model], DataFrame)
    #         rename!(adm, (n => split(n)[1] for n in names(adm))...)
    #         @test adm == ustrip(extract_sac(HEADER1, COLUMNS1) |> age_depth_model)
    #     end
    # end

    # @testset "Waterdepth signs" begin
    #     BS92_TEST_INPUT = BS92.Input(
    #         tag = "single pixel model",
    #         box = Box{Periodic{2}}(grid_size=(1, 1), phys_scale=600.0u"m"),
    #         time = TimeProperties(
    #           Δt = 10.0u"yr",
    #           steps = 8000),
    #         output = Dict(:full => OutputSpec((1, 1), 80)),
    #         sea_level = t -> 10.0u"m" * sin(2π * t / 20u"kyr"),
    #         initial_topography = (_, _) -> - 50.0u"m",
    #         subsidence_rate = 0.001u"m/yr",
    #         insolation = 400.0u"W/m^2",
    #         facies = [BS92.Facies(
    #           maximum_growth_rate = 0.005u"m/yr",
    #           saturation_intensity = 50.0u"W/m^2",
    #           extinction_coefficient = 0.05u"m^-1"
    #         )])

    #     mktempdir() do path
    #         run_model(Model{BS92}, BS92_TEST_INPUT, joinpath(path, "run.h5"))
    #         header, data = read_column(joinpath(path, "run.h5"), :full)
    #         wd = extract_wd(header, data, 1)
    #         sac = extract_sac(header, data, 1)
    #         submerged = wd.wd1 .> -1.0u"m"
    #         growing = (sac.sac_1[2:end] .- sac.sac_1[1:end-1]) .> 0.5u"m"
    #         @test all(growing .&& (submerged[1:end-1] .|| submerged[2:end]) .|| .!growing)
    #     end
    # end
end
# ~/~ end
