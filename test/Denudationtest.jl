using Test 
using Unitful

@testset "DenudationTST" begin
    import CarboKitten.Denudation.CarbDissolution: dissolution
    import CarboKitten.Denudation.EmpericalDenudation: emperical_denudation, slope_kernel
    import CarboKitten.Denudation.PhysicalErosion: physical_erosion, mass_erosion, total_mass_redistribution
    using CarboKitten.Stencil: Periodic, Reflected, stencil
    using CarboKitten.Config: Box, Vectors, TimeProperties
    using CarboKitten.Burgess2013.CA: step_ca, run_ca
    using CarboKitten.Burgess2013.Config: Facies
    using CarboKitten.InputConfig: Input, DenudationType
    using CarboKitten.Denudation: denudation, calculate_redistribution, Dissolution, NoDenudation, PhysicalErosionParam, EmpericalDenudationParam
    

    DENUDATION_LOW_T = Dissolution(273.0,1000.0,10^(-1.5),2e-3)
    DENUDATION_HIGH_T = Dissolution(303.0,1000.0,10^(-1.5),2e-3)
    DENUDATION_LOW_P = EmpericalDenudationParam(800.0)
    DENUDATION_HIGH_P = EmpericalDenudationParam(1000.0)
    DENUDATION_PHYS = PhysicalErosionParam(0.23)
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

    box = Box{Periodic{2}}(grid_size=(5, 5), phys_scale=100.0u"m")
    n_facies = length(MODEL1)
    ca_init = rand(0:n_facies, box.grid_size...)
    ca_init = [ 0  0  1  3  3
                0  1  3  2  1
                2  0  1  0  1
                1  3  3  3  0
                1  3  2  3  2]
    denudation_mass_LOW_T = zeros(typeof(0.0u"m/kyr"),box.grid_size...)
    denudation_mass_HIGH_T = zeros(typeof(0.0u"m/kyr"),box.grid_size...)
    denudation_mass_LOW_P = zeros(typeof(0.0u"m/kyr"),box.grid_size...)
    denudation_mass_HIGH_P = zeros(typeof(0.0u"m/kyr"),box.grid_size...)
    denudation_mass_phys = zeros(typeof(0.0u"m/kyr"),box.grid_size...)
    denudation_mass_phys_flat = zeros(typeof(0.0u"m/kyr"),box.grid_size...)
    redistribution_mass = zeros(typeof(0.0u"m/kyr"),box.grid_size...)
    water_depth = rand(box.grid_size...)
    water_depth = [ 0.989943  0.48076   0.518983  0.997996   0.895681
                    0.872733  0.208779  0.882917  0.550494   0.674066
                    0.57987   0.619433  0.769506  0.593786   0.856186
                    0.407728  0.469545  0.896348  0.473817   0.797112
                    0.610194  0.921632  0.322729  0.0103646  0.691191]
    water_depth_flat = 0.5 .* ones(box.grid_size...)
    inf_map = ones(box.grid_size...)
    slope = rand(Float64, box.grid_size...)
    slopefn = stencil(Float64, Periodic{2}, (3, 3), slope_kernel)
    slopefn(water_depth, slope, box.phys_scale ./u"m")
    slope_flat = zeros(box.grid_size...)
    slopefn(water_depth_flat, slope_flat, box.phys_scale ./u"m")
    for idx in CartesianIndices(ca_init)
        f = ca_init[idx]
        if f == 0
            continue
        end
    (denudation_mass_LOW_T[idx]) = denudation(box, DENUDATION_LOW_T, water_depth[idx], slope[idx],MODEL1[f])
    (denudation_mass_HIGH_T[idx]) = denudation(box, DENUDATION_HIGH_T, water_depth[idx], slope[idx],MODEL1[f])
    end


    for idx in CartesianIndices(ca_init)
        f = ca_init[idx]
        if f == 0
            continue
        end
    (denudation_mass_LOW_P[idx]) = denudation(box, DENUDATION_LOW_P, water_depth[idx], slope[idx],MODEL1[f])
    (denudation_mass_HIGH_P[idx]) = denudation(box, DENUDATION_HIGH_P, water_depth[idx], slope[idx],MODEL1[f])
    end

    for idx in CartesianIndices(ca_init)
        f = ca_init[idx]
        if f == 0
            f = 1
        end
    (denudation_mass_phys[idx]) = denudation(box, DENUDATION_PHYS, water_depth[idx], slope[idx],MODEL1[f])
    (denudation_mass_phys_flat[idx]) = denudation(box, DENUDATION_PHYS, water_depth_flat[idx], slope_flat[idx],MODEL1[f])
    inf_map[idx] = MODEL1[f].infiltration_coefficient
    end
    (redistribution_mass) = calculate_redistribution(box,DENUDATION_PHYS,water_depth,slope,inf_map)

    println("Denudation mass: ", denudation_mass_phys ./u"m/kyr")
    println("Redistribution mass: ", redistribution_mass ./u"m/kyr")
    println("Denudation mass, Low T: ", denudation_mass_LOW_T ./u"m/kyr")
    println("Denudation mass, Low P: ", denudation_mass_LOW_P ./u"m/kyr")
    @test sum(denudation_mass_LOW_T) < sum(denudation_mass_HIGH_T) 
    @test sum(denudation_mass_LOW_P) < sum(denudation_mass_HIGH_P)
    @test sum(denudation_mass_phys) > sum(denudation_mass_phys_flat)
    @test sum(denudation_mass_phys) ≈ sum(redistribution_mass) 

    #regression_test
    @test 1651*0.95 < abs.(sum(denudation_mass_LOW_T)) ./u"m/kyr" < 1651 *1.05
    @test 479.5*0.95 < abs.(sum(denudation_mass_LOW_P)) ./u"m/kyr" < 479.5*1.05
    @test 0.95* 0.7115 < sum(denudation_mass_phys) ./u"m/kyr" < 0.7115*1.05

end