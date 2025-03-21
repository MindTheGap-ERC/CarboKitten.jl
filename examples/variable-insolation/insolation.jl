using CarboKitten

module VariableInsolation

using CarboKitten
using Unitful
using DelimitedFiles: readdlm

const PATH = "data/output"

const TAG = "insolation-example"

function insolation()
    dir = @__DIR__
    filename = joinpath(dir, "insolation.txt")
    insol = readdlm(filename, '\t', header=false) 
    return insol[2:end] .* u"W/m^2"
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

const PERIOD = 0.2u"Myr"
const AMPLITUDE = 0.5u"m"

const INPUT = ALCAP.Input(
    tag="$TAG",
    box=Box{Coast}(grid_size=(100, 50), phys_scale=170.0u"m"),
    time=TimeProperties(
        Δt=200.0u"yr",
        steps=1000,
        write_interval=1),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=t -> AMPLITUDE * sin(2π * t / PERIOD),
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    insolation=insolation(),
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    facies=FACIES)

end

run_model(Model{ALCAP}, VariableInsolation.INPUT, "data/output/varins.h5")

