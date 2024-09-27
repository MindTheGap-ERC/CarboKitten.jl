# ~/~ begin <<docs/src/model-alcap.md#examples/alcaps/alternative.jl>>[init]
#| requires: src/Model/ALCAPS.jl
#| creates: data/alcaps2.h5

using Unitful
using CarboKitten.BoundaryTrait: Shelf
using CarboKitten.Config: Box, TimeProperties
using CarboKitten.Model.ALCAPS: Facies, Input, main
using CarboKitten.Export: data_export, CSV

const m = u"m"
const Myr = u"Myr"

const PATH = "data"
const TAG = "alcaps2"

const FACIES = [
    Facies(viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=10000u"m"),
    Facies(viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=5000u"m"),
    Facies(viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=7000u"m")
]

const PERIOD = 0.2Myr
const AMPLITUDE = 4.0m

const INPUT = Input(
    tag="ALCAPS alternative",
    box=Box{Shelf}(grid_size=(100, 50), phys_scale=150.0m),
    time=TimeProperties(
        Δt=0.0002Myr,
        steps=5000,
        write_interval=1),
    ca_interval=1, bedrock_elevation=(x, y) -> -x / 300.0,
    sea_level=t -> AMPLITUDE * sin(2π * t / PERIOD),
    subsidence_rate=50.0m / Myr,
    disintegration_rate=500.0m / Myr,
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5m,
    facies=FACIES)

main(INPUT, "$(PATH)/alcaps2.h5")

data_export(
    CSV(tuple.(10:20:70, 25),
      :sediment_accumulation_curve => "$(PATH)/$(TAG)_sac.csv",
      :age_depth_model => "$(PATH)/$(TAG)_adm.csv",
      :stratigraphic_column => "$(PATH)/$(TAG)_sc.csv",
      :metadata => "$(PATH)/$(TAG).toml"),
    "$(PATH)/alcaps2.h5")
# ~/~ end
