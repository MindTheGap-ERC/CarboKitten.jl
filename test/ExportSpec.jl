# ~/~ begin <<docs/src/data-export.md#test/ExportSpec.jl>>[init]
using CarboKitten
using CarboKitten.Export: Axes, Header, DataVolume, data_export, CSVExportTrait,
    age_depth_model, extract_sac, extract_sc, CSV, read_data, extract_sac, extract_wd,
    read_column
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
    grid_size=(3, 1),
    n_facies=1,
    initial_topography=zeros(typeof(1.0u"m"), 3, 3),
    sea_level=zeros(typeof(1.0u"m"), 11),
    subsidence_rate=10u"m/Myr",
    data_sets=Dict())

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
    deposition=PRODUCTION1,
    sediment_thickness=ELEVATION1)

const GRID_LOCATIONS1 = [(1, 1), (2, 1), (3, 1)]

const COLUMNS1 = [DATA1[loc...] for loc in GRID_LOCATIONS1]
# ~/~ end

@testset "CarboKitten.Export" begin
    sac = data_export(CSVExportTrait{:sediment_accumulation_curve}, HEADER1, COLUMNS1)
    adm = data_export(CSVExportTrait{:age_depth_model}, HEADER1, COLUMNS1)
    sc = data_export(CSVExportTrait{:stratigraphic_column}, HEADER1, COLUMNS1)

    # ~/~ begin <<docs/src/data-export.md#export-test>>[init]
    @testset "Hither and Dither" begin
        @test sac.sac_1 ≈ ELEVATION1[1, 1, :]
        @test sac.sac_2 ≈ ELEVATION1[2, 1, :]
        @test sac.sac_3 ≈ ELEVATION1[3, 1, :]
    end
    # ~/~ end
    # ~/~ begin <<docs/src/data-export.md#export-test>>[1]
    @testset "ADM Monotonicity" begin
        @test sac.sac_1 == adm.adm_1
        @test sac.sac_3 == adm.adm_3
        @test sac.sac_2 != adm.adm_2

        @test all(adm.adm_2[2:end] .- adm.adm_2[1:end-1] .>= 0.0u"m")
    end
    # ~/~ end
    # ~/~ begin <<docs/src/data-export.md#export-test>>[2]
    @testset "SC sum equals ADM" begin
        @test [0.0u"m"; cumsum(sc.sc_1_f1)] ≈ adm.adm_1
        @test [0.0u"m"; cumsum(sc.sc_2_f1)] ≈ adm.adm_2
        @test [0.0u"m"; cumsum(sc.sc_3_f1)] ≈ adm.adm_3
    end
    # ~/~ end
    @testset "Write to folder" begin
        using DataFrames: select

        mktempdir() do path
            spec = CSV(
                :sediment_accumulation_curve => joinpath(path, "sac.csv"),
                :age_depth_model => joinpath(path, "adm.csv"),
                :stratigraphic_column => joinpath(path, "sc.csv"),
                :water_depth => joinpath(path, "wd.csv"),
                :metadata => joinpath(path, "metadata.toml"))
            data_export(spec, HEADER1, COLUMNS1)
            for f in values(spec.output_files)
                @test isfile(f)
            end

            metadata = TOML.parsefile(spec.output_files[:metadata])
            @test IdDict(Symbol(k) => v for (k, v) in metadata["files"]) == spec.output_files
            @test length(metadata["locations"]) == 3

            adm_tab = read_csv(spec.output_files[:age_depth_model], DataFrame)
            rename!(adm_tab, (n => split(n)[1] for n in names(adm_tab))...)
            @test select(adm_tab, ["adm_$(i)" for i in 1:3]) == 
                select(ustrip(adm), ["adm_$(i)" for i in 1:3])
        end
    end

    @testset "Water depth signs" begin
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