using Test 
using Unitful

@testset "DenudationTST" begin
    import CarboKitten.Denudation.CarbDissolution: dissolution
    import CarboKitten.Denudation.EmpericalDenudation: emperical_denudation, slope_kernel
    import CarboKitten.Denudation.PhysicalErosion: physical_erosion, mass_erosion, total_mass_redistribution
    using CarboKitten.Stencil: Periodic, stencil
    using CarboKitten.Config: Box, Vectors, TimeProperties
    using CarboKitten.Burgess2013.CA: step_ca, run_ca
    using CarboKitten.Burgess2013.Config: Facies
    using CarboKitten.InputConfig: Input, DenudationType
    using CarboKitten.Denudation: denudation, Dissolution, NoDenudation, PhysicalErosionParam, EmpericalDenudationParam
    

    DENUDATION_LOW_T = Dissolution(273.0,100.0,10^(-1.5),2e-3)
    DENUDATION_HIGH_T = Dissolution(303.0,1000.0,10^(-1.5),2e-3)

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

    box = Box{Periodic{2}}(grid_size=(50, 50), phys_scale=100.0u"m")
    n_facies = length(MODEL1)
    ca_init = rand(0:n_facies, box.grid_size...)
    denudation_mass_LOW_T = zeros(typeof(0.0u"m/kyr"),box.grid_size...)
    denudation_mass_HIGH_T = zeros(typeof(0.0u"m/kyr"),box.grid_size...)
    water_depth = rand(box.grid_size...)

    slope = zeros(Float64, box.grid_size...)
    slopefn = stencil(Float64, Periodic{2}, (3, 3), slope_kernel)
    slopefn(water_depth, slope, box.phys_scale ./u"m")

    for idx in CartesianIndices(ca_init)
        f = ca_init[idx]
        if f == 0
            continue
        end
    (denudation_mass_LOW_T[idx]) = denudation(box, DENUDATION_LOW_T, water_depth[idx], slope[idx],MODEL1[f])
    (denudation_mass_HIGH_T[idx]) = denudation(box, DENUDATION_HIGH_T, water_depth[idx], slope[idx],MODEL1[f])
    end

    @test sum(denudation_mass_LOW_T) < sum(denudation_mass_HIGH_T) 


end