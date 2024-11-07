# ~/~ begin <<docs/src/model-alcap.md#examples/model/alcap/run.jl>>[init]
#| requires: src/Model/ALCAP2.jl
#| creates: data/output/alcap2.h5

module Script

using Unitful
using CarboKitten.Components
using CarboKitten.Components.Common
using CarboKitten.Model: ALCAP2 as ALCAP
using CarboKitten.Export: data_export, CSV
using DataFrames
using XLSX
using Interpolations

const m = u"m"
const Myr = u"Myr"

const PATH = "data/output"
const TAG = "alcap2"
const PERIOD = 0.2Myr
const AMPLITUDE = 4.0m
const filepath = "data/input/La04.xlsx"

function sealevel_curve(filepath)
    table = DataFrame(XLSX.readtable(filepath, "T1")) 
    time = table.Yr ./ 1e6 .* u"Myr"
    sl = table.SL .* 1.0u"m"

    data = DataFrame(time = time, sea_level = sl)

    x  = linear_interpolation(data.time, data.sea_level)
    return x
end 



const FACIES = [
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=10000u"m"),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=5000u"m"),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=7000u"m")
]

const INPUT = ALCAP.Input(
    tag="$TAG",
    box=Box{Shelf}(grid_size=(100, 50), phys_scale=150.0m),
    time=TimeProperties(
        t0 = -2.0u"Myr",
        Δt=0.001Myr,
        steps=1500,
        write_interval=1),
    ca_interval=1,
    bedrock_elevation=(x, y) -> -x / 300.0,
    sea_level= sealevel_curve(filepath),
    subsidence_rate=50.0m / Myr,
    disintegration_rate=500.0m / Myr,
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5m,
    facies=FACIES)

function main()
    H5Writer.run(Model{ALCAP}, INPUT, "$(PATH)/$(TAG).h5")

    data_export(
        CSV(tuple.(10:20:70, 25),
          :sediment_accumulation_curve => "$(PATH)/$(TAG)_sac.csv",
          :age_depth_model => "$(PATH)/$(TAG)_adm.csv",
          :stratigraphic_column => "$(PATH)/$(TAG)_sc.csv",
          :metadata => "$(PATH)/$(TAG).toml"),
        "$(PATH)/$(TAG).h5")
end

end

Script.main()
# ~/~ end