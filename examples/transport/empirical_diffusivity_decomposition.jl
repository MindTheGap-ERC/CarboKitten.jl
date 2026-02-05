# ~/~ begin <<docs/src/active-layer-transport.md#examples/transport/empirical_diffusivity_decomposition.jl>>[init]
include("Kaufman_diffusivity.jl")
using .KaufmanCommon

using CarboKitten
using CarboKitten.Models.ALCAP
using Unitful
using GeometryBasics

# watch this use of a closure! Gives an end error somewhere but it runs
function make_input_decomposition(;
        grid_size = (100, 1),
        phys_scale = 150.0u"m",
        Δt = 0.0002u"Myr",
        steps = 2500, 
        tag = "kaufman_decomposition")

    D_wave_surface = KaufmanCommon.KAUFMAN_C0 - KaufmanCommon.BACKGROUND_DIFFUSION
    T_wave = 8.0u"s"
    omega = (2π / T_wave) |> ustrip |> x -> x * u"1/yr"

    # Back-calculate surface orbital velocity 
    # u0 = sqrt(D_wave_surface * omega / (u"m/yr"))
    # this cannot be applied because of unit mismatch - why does CK use diffusivity units m/yr? fixme
    # Meanwhile we set this to some constant number that will hopefully prove realistic
    u0 = 100.0u"m/yr"
    # Wave decay constant (half of Kaufman's because velocity not diffusivity) 
    k_wave = KaufmanCommon.KAUFMAN_C1 / 2

    wave_velocity_fn = function(water_depth)
        u = u0 * exp(-k_wave * water_depth)
        # Onshore transport
        v = Vec2(-u, 0.0u"m/yr") 
        # Shear from exponential decay 
        s = Vec2(k_wave * u, 0.0u"1/yr")  
        return (v, s)
    end

    facies = [
        ALCAP.Facies(
            viability_range = (4, 10),
            activation_range = (6, 10),
            maximum_growth_rate = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2",
            diffusion_coefficient = KaufmanCommon.BACKGROUND_DIFFUSION/2,
            wave_velocity = wave_velocity_fn,
            name = "producer"),
        ALCAP.Facies(
            active = false,
            diffusion_coefficient = KaufmanCommon.BACKGROUND_DIFFUSION/2,
            wave_velocity = wave_velocity_fn,
            name = "transported")
    ]

    PERIOD = 0.2u"Myr"
    AMPLITUDE = 4.0u"m"

    return ALCAP.Input(
        tag = tag,
        box = CarboKitten.Box{Coast}(grid_size=grid_size, phys_scale=phys_scale),
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
    input = make_input_decomposition()
    output = "data/output/kaufman_decomposition.h5"
    run_model(Model{ALCAP}, input, output)
end

main()
end
# ~/~ end
