# ~/~ begin <<docs/src/cases/tabular-sea-level.md#examples/tabular-sea-level/run.jl>>[init]

module TabularSeaLevel

using CarboKitten

using DelimitedFiles: readdlm
using DataFrames
using CarboKitten.DataSets: artifact_dir
using Interpolations
using CategoricalArrays

# ~/~ begin <<docs/src/cases/tabular-sea-level.md#tabular-sea-level>>[init]
function miller_2020()
    dir = artifact_dir()
    filename = joinpath(dir, "Miller2020", "Cenozoic_sea_level_reconstruction.tab")

    data, header = readdlm(filename, '\t', header=true)
    return DataFrame(
        time=-data[:,4] * u"kyr",
        sealevel=data[:,7] * u"m",
        refkey=categorical(data[:,2]),
        reference=categorical(data[:,3]))
end
# ~/~ end
# ~/~ begin <<docs/src/cases/tabular-sea-level.md#tabular-sea-level>>[1]
function sea_level()
    df = miller_2020()
    lisiecki_df = df[df.refkey .== "846 Lisiecki", :]
    sort!(lisiecki_df, [:time])

    return linear_interpolation(
        lisiecki_df.time,
        lisiecki_df.sealevel)
end
# ~/~ end
# ~/~ begin <<docs/src/cases/tabular-sea-level.md#tabular-sea-level>>[2]
const TIME_PROPERTIES = TimeProperties(
    t0 = -2.0u"Myr",
    Δt = 200.0u"yr",
    steps = 5000
)
# ~/~ end
# ~/~ begin <<docs/src/cases/tabular-sea-level.md#tabular-sea-level>>[3]
const PATH = "data/output"
const TAG = "lisiecki-sea-level"

const FACIES = [
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=50.0u"m/yr"),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=25.0u"m/yr"),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=35.0u"m/yr")
]

const INPUT = ALCAP.Input(
    tag="$TAG",
    box=CarboKitten.Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m"),
    time=TIME_PROPERTIES,
    ca_interval=1,
    initial_topography=(x, y) -> -x / 200.0,
    sea_level=sea_level(),
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    facies=FACIES)

function main()
    run_model(Model{ALCAP}, INPUT, "$(PATH)/$(TAG).h5")
end
# ~/~ end

end

TabularSeaLevel.main()
# ~/~ end
