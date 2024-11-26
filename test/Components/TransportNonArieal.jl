module TransportNonArieal
using Test
using CarboKitten.Components.Common
using CarboKitten.Components.WaterDepth: water_depth
using CarboKitten.Components.ActiveLayer: disintegration
using CarboKitten.Models.ALCAP: Facies, Input, initial_state
using Unitful

@testset "Components/TransportNonArieal" begin
    let FACIES = [
        Facies(
            viability_range=(4, 10),
            activation_range=(6, 10),
            maximum_growth_rate=500u"m/Myr",
            extinction_coefficient=0.8u"m^-1",
            saturation_intensity=60u"W/m^2",
            diffusion_coefficient=10000u"m"),
        Facies(
            viability_range=(4, 10),
            activation_range=(6, 10),
            maximum_growth_rate=400u"m/Myr",
            extinction_coefficient=0.1u"m^-1",
            saturation_intensity=60u"W/m^2",
            diffusion_coefficient=5000u"m"),
        Facies(
            viability_range=(4, 10),
            activation_range=(6, 10),
            maximum_growth_rate=100u"m/Myr",
            extinction_coefficient=0.005u"m^-1",
            saturation_intensity=60u"W/m^2",
            diffusion_coefficient=7000u"m")
    ]
    input = Input(
        box=Box{Coast}(grid_size=(5, 5), phys_scale=150.0u"m"),
        time=TimeProperties(
            Δt=0.0002u"Myr",
            steps=5000,
            write_interval=1),
        ca_interval=1,
        initial_topography=(x, y) -> x / 300.0,
        sea_level= t -> 0.0u"m",
        subsidence_rate=0.0u"m/Myr",
        disintegration_rate=50.0u"m/Myr",
        insolation=400.0u"W/m^2",
        sediment_buffer_size=5,
        depositional_resolution=0.5u"m",
        facies=FACIES)
    state = initial_state(input)
    state.sediment_height = 20 .* rand(5,5) *u"m"
    state.sediment_buffer = rand(Float64, input.sediment_buffer_size, 3 , input.box.grid_size...)
    rs = disintegration(input)
    size(rs(state))
    @test sum(rs(state)) == 0.0*u"m"
    end
end
end
