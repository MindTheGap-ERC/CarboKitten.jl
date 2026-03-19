# ~/~ begin <<docs/src/active-layer-transport.md#examples/model/diffusivity/alcap_diffusivity.jl>>[init]
module Diffusivity_example

using Unitful
using CarboKitten
using CarboKitten.Export: read_slice, data_export, CSV

const PATH = "data/output"

const TAG = "diffusivity-example"

cost min_diffusivity = 2.5u"m/yr"

const FACIES = [
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        transport_coefficient=4*min_diffusivity,
        name="euphotic"),
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        transport_coefficient=min_diffusivity,
        name="oligophotic"),
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        transport_coefficient=2*min_diffusivity,
        name="aphotic"),

    ALCAP.Facies(
        active=false,
        transport_coefficient=10.0u"m/yr",
        name="euphotic transported"),
    ALCAP.Facies(
        active=false,
        transport_coefficient=1.0u"m/yr",
        name="oligophotic transported"),
    ALCAP.Facies(
        active=false,
        transport_coefficient=5.0u"m/yr",
        name="aphotic transported")
]

const PERIOD = 0.2u"Myr"
const AMPLITUDE = 4.0u"m"

const INPUT = ALCAP.Input(
    tag="$TAG",
    box=Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m"),
    time=TimeProperties(
        Δt=0.0002u"Myr",
        steps=5000),
    output=Dict(
        :topography => OutputSpec(slice=(:,:), write_interval=10),
        :profile => OutputSpec(slice=(:, 25), write_interval=1)),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=t -> AMPLITUDE * sin(2π * t / PERIOD),
    subsidence_rate=50.0u"m/Myr",
    cementation_time = 5000u"yr",
    disintegration_rate=5.0u"m/Myr",
    disintegration_transfer = f -> stack((0.0.*f[1,:,:], 0.0.*f[2,:,:], 0.0.*f[3,:,:],
                                      f[1,:,:].+f[4,:,:], f[2,:,:].+f[5,:,:], f[3,:,:].+f[6,:,:]), dims=1),
    insolation=400.0u"W/m^2",
    sediment_buffer_size=2,
    depositional_resolution=0.5u"m",
    facies=FACIES)

function main()
    run_model(Model{ALCAP}, INPUT, "$(PATH)/$(TAG).h5")
end

end

Diffusivity_example.main()

# ~/~ end
