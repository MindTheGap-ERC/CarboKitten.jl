# ~/~ begin <<docs/src/active-layer-transport.md#examples/transport/empirical_diffusivity_substitution.jl>>[init]
include("Kaufman_diffusivity.jl")
using .KaufmanCommon

using CarboKitten
using CarboKitten.Models.ALCAP
using Unitful
using GeometryBasics

function make_input_substitution(;
        grid_size = (100, 1),
        phys_scale = 150.0u"m",
        Δt = 0.0002u"Myr",
        steps = 2500,
        tag = "kaufman_substitution")

    avg_depth = 5.0u"m"
    D_eff = KaufmanCommon.kaufman_diffusivity(avg_depth)

    facies = [
        ALCAP.Facies(
            viability_range = (4, 10),
            activation_range = (6, 10),
            maximum_growth_rate = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2",
            diffusion_coefficient = D_eff,
            wave_velocity = _ -> (Vec2(0.0, 0.0)u"m/yr", Vec2(0.0, 0.0)u"1/yr"),
            name = "producer"),
        ALCAP.Facies(
            active = false,
            diffusion_coefficient = D_eff,
            wave_velocity = _ -> (Vec2(0.0, 0.0)u"m/yr", Vec2(0.0, 0.0)u"1/yr"),
            name = "transported")
    ]

    PERIOD = 0.2u"Myr"
    AMPLITUDE = 4.0u"m"

    return ALCAP.Input(
        tag = tag,
        box = Box{Coast}(grid_size=grid_size, phys_scale=phys_scale),
        time = TimeProperties(Δt=Δt, steps=steps),
        output = Dict(
            :profile => OutputSpec(slice=(:, 1), write_interval=10)),
        ca_interval = 1,
        initial_topography = (x, y) -> -x / 300.0,
        sea_level = t -> AMPLITUDE * sin(2π * t / PERIOD),
        subsidence_rate = 50.0u"m/Myr",
        disintegration_rate = 50.0u"m/Myr",
        disintegration_transfer = p -> [0.0u"m", p[1] + p[2]],
        insolation = 400.0u"W/m^2",
        sediment_buffer_size = 50,
        depositional_resolution = 0.5u"m",
        facies = facies)
end

function main()
    KaufmanCommon.print_diffusivity_profile()
    input = make_input_substitution()
    output = "data/output/kaufman_substitution.h5"
    run_model(Model{ALCAP}, input, output)
end

main()
# ~/~ end
