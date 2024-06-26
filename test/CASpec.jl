# ~/~ begin <<docs/src/carbocat-ca.md#test/CASpec.jl>>[init]
using Test 
@testset "CA" begin
    using CarboKitten.BoundaryTrait: Periodic
    using CarboKitten.Config: Box
    using CarboKitten.Burgess2013.CA: step_ca, run_ca
    using CarboKitten.Burgess2013.Config: Facies
    using Unitful

    mutable struct State
        ca::Matrix{Int}
        ca_priority::Vector{Int}
    end

    MODEL1 = [
        Facies(viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate = 500u"m/Myr",
        extinction_coefficient = 0.8u"m^-1",
        saturation_intensity = 60u"W/m^2",
        reactive_surface = 1000,
        mass_density = 2730,
        infiltration_coefficient= 0.5),
    
        Facies(viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate = 400u"m/Myr",
        extinction_coefficient = 0.1u"m^-1",
        saturation_intensity = 60u"W/m^2",
        reactive_surface = 1000,
        mass_density = 2730,
        infiltration_coefficient= 0.5),
    
        Facies(viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate = 100u"m/Myr",
        extinction_coefficient = 0.005u"m^-1",
        saturation_intensity = 60u"W/m^2",
        reactive_surface = 1000,
        mass_density = 2730,
        infiltration_coefficient= 0.5)
    ]

    n_facies = length(MODEL1)
    box = Box{Periodic{2}}(grid_size=(50, 50), phys_scale=100.0u"m")
    ca_init = rand(0:n_facies, box.grid_size...)
    ca_channel = run_ca(box, MODEL1, copy(ca_init), n_facies)
    item1, _ = Iterators.peel(Iterators.drop(ca_channel, 20))

    ca_step = step_ca(box, MODEL1)
    state = State(copy(ca_init), 1:n_facies)
    for _ in 1:20
        ca_step(state)
    end

    @test item1 == state.ca
end
# ~/~ end
