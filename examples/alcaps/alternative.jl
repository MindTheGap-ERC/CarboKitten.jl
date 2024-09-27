# ~/~ begin <<docs/src/model-alcap.md#examples/alcaps/alternative.jl>>[init]
#| requires: src/Model/ALCAPS.jl
#| creates: data/alcaps2.h5

using Unitful
using CarboKitten.BoundaryTrait: Shelf
using CarboKitten.Config: Box, TimeProperties
using CarboKitten.Model.ALCAPS: Facies, Input, main

const m = u"m"
const Myr = u"Myr"

const FACIES = [
    Facies(viability_range = (4, 10),
           activation_range = (6, 10),
           maximum_growth_rate = 500u"m/Myr",
           extinction_coefficient = 0.8u"m^-1",
           saturation_intensity = 60u"W/m^2",
           diffusion_coefficient = 10000u"m"),

    Facies(viability_range = (4, 10),
           activation_range = (6, 10),
           maximum_growth_rate = 400u"m/Myr",
           extinction_coefficient = 0.1u"m^-1",
           saturation_intensity = 60u"W/m^2",
           diffusion_coefficient = 5000u"m"),

    Facies(viability_range = (4, 10),
           activation_range = (6, 10),
           maximum_growth_rate = 100u"m/Myr",
           extinction_coefficient = 0.005u"m^-1",
           saturation_intensity = 60u"W/m^2",
           diffusion_coefficient = 7000u"m")
]

const PERIOD = 0.2Myr
const AMPLITUDE = 4.0m

const INPUT = Input(
    box = Box{Shelf}(grid_size=(100, 50), phys_scale=150.0m),
    time = TimeProperties(
        Δt=0.0002Myr,
        steps=5000,
        write_interval=1),
    ca_interval           = 1,

    bedrock_elevation        = (x, y) -> -x / 300.0 ,
    sea_level                = t -> AMPLITUDE * sin(2π * t / PERIOD),
    subsidence_rate          = 50.0m/Myr,
    disintegration_rate      = 500.0m/Myr,
    insolation               = 400.0u"W/m^2",
    sediment_buffer_size     = 50,
    depositional_resolution  = 0.5m,
    facies                   = FACIES)

main(INPUT, "data/alcaps2.h5")
# ~/~ end