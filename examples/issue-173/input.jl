module Script

using Unitful
using CarboKitten

const PATH = "data"

const TAG = "carbonate_stratpal_1"

const FACIES = [
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=50.0u"m/yr"),
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=25.0u"m/yr"),
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=12.5u"m/yr")
]

const PERIOD1 = 2.0u"Myr"
const AMPLITUDE1 = 15.0u"m"
const PERIOD2 = 0.2u"Myr"
const AMPLITUDE2 = 2.5u"m"

const INPUT = ALCAP.Input(
    diagnostics=true,
    tag="$TAG",
    box=Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m"),
    time=TimeProperties(
        Δt=100u"yr",
        steps=40000),
    output=Dict(
        :topgraphy => OutputSpec(write_interval=1000),
        :profile => OutputSpec(slice=(:, 25), write_interval=10)),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=t -> AMPLITUDE1 * sin(2π * t / PERIOD1) + AMPLITUDE2 * sin(2π * t / PERIOD2),
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    cementation_time=100.0u"yr",
    facies=FACIES)


function main()
    run_model(Model{ALCAP}, INPUT, "$(PATH)/$(TAG).h5")

end

end

Script.main()
