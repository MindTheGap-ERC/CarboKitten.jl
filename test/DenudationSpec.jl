module DenudationSpec
using Test
using Unitful

using CarboKitten
using CarboKitten.Denudation.DissolutionMod: dissolution, dominant_facies
using CarboKitten.SedimentStack: push_sediment!, peek_sediment
using CarboKitten.Components.Common
using CarboKitten.Components: Denudation as D


#import CarboKitten.Denudation.PhysicalErosionMod: physical_erosion, mass_erosion, total_mass_redistribution
using CarboKitten.Stencil: Periodic, Reflected, stencil
using CarboKitten.Config: Vectors, TimeProperties
using CarboKitten.Boxes: Box
using CarboKitten.Models: WithDenudation as WD
using CarboKitten.Denudation.EmpiricalDenudationMod: slope_kernel
using CarboKitten.Denudation: denudation, redistribution, Dissolution, NoDenudation, PhysicalErosion, EmpiricalDenudation


@testset "DenudationTST" begin
 
    DENUDATION_HIGH_CO2 = Dissolution(temp = 293.0u"K",precip = 1.0u"m/yr", pco2 = 10^(-1.5)*u"atm",reactionrate = 2e-3u"m/yr")
    DENUDATION_LOW_CO2 = Dissolution(temp = 293.0u"K",precip = 1.0u"m/yr", pco2 = 10^(-2.5)*u"atm",reactionrate = 2e-3u"m/yr")
    DENUDATION_LOW_P = EmpiricalDenudation(precip = 0.8u"m/yr")
    DENUDATION_HIGH_P = EmpiricalDenudation(precip = 1.0u"m/yr")
    DENUDATION_PHYS = PhysicalErosion()
    
    MODEL1 = [
        WD.Facies(viability_range = (4, 10),
        activation_range = (6, 10),
        production = BenthicProduction(
            maximum_growth_rate = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        transport_coefficient=50u"m/yr",
        reactive_surface = 1000u"m^2/m^3",
        mass_density = 2730u"kg/m^3",
        infiltration_coefficient= 0.5,
        erodibility = 0.23u"m/yr"),

        WD.Facies(viability_range = (4, 10),
        activation_range = (6, 10),
        production = BenthicProduction(
            maximum_growth_rate = 400u"m/Myr",
            extinction_coefficient = 0.1u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        transport_coefficient= 50u"m/yr",
        reactive_surface = 1000u"m^2/m^3",
        mass_density = 2730u"kg/m^3",
        infiltration_coefficient= 0.5,
        erodibility = 0.23u"m/yr"),

        WD.Facies(viability_range = (4, 10),
        activation_range = (6, 10),
        production = BenthicProduction(
            maximum_growth_rate = 100u"m/Myr",
            extinction_coefficient = 0.005u"m^-1",
            saturation_intensity = 60u"W/m^2"),
        transport_coefficient= 50u"m/yr",
        reactive_surface = 1000u"m^2/m^3",
        mass_density = 2730u"kg/m^3",
        infiltration_coefficient= 0.5,
        erodibility = 0.23u"m/yr")
    ]

    box = Box{Periodic{2}}(grid_size=(5, 5), phys_scale=1.0u"km")
    n_facies = length(MODEL1)
    ca_init = [ 0  0  1  3  3
                0  1  3  2  1
                2  0  1  0  1
                1  3  3  3  0
                1  3  2  3  2]


    STATE1 = WD.State(
        step=0,
        bathymetry=zeros(Height, box.grid_size...),
        sediment_thickness=zeros(Height, box.grid_size...),
        sediment_buffer=zeros(Float64, 10, n_facies, box.grid_size...),
        active_layer=zeros(Amount, n_facies, box.grid_size...),
        ca=ca_init, ca_priority=[1,2,3]
    )

    denudation_mass_HIGH_CO2 = zeros(typeof(0.0u"m/kyr"),n_facies,box.grid_size...)
    denudation_mass_LOW_CO2 = zeros(typeof(0.0u"m/kyr"),n_facies,box.grid_size...)
    denudation_mass_LOW_P = zeros(typeof(0.0u"m/kyr"),n_facies,box.grid_size...)
    denudation_mass_HIGH_P = zeros(typeof(0.0u"m/kyr"),n_facies,box.grid_size...)
    denudation_mass_phys = zeros(typeof(0.0u"m/kyr"),n_facies,box.grid_size...)
    denudation_mass_phys_flat = zeros(typeof(0.0u"m/kyr"),n_facies,box.grid_size...)
    redistribution_mass = zeros(typeof(0.0u"m"),n_facies,box.grid_size...)

    water_depth = -100 .* [ 0.989943  0.48076   0.518983  0.997996   0.895681;
                    0.872733  0.208779  0.882917  0.550494   0.674066;
                    0.57987   0.619433  0.769506  0.593786   0.856186;
                    0.407728  0.469545  0.896348  0.473817   0.797112;
                    0.610194  0.921632  0.322729  0.0103646  0.691191]

    water_depth_flat = -0.5 .* ones(box.grid_size...)
    inf_map = ones(box.grid_size...)
    slope = rand(Float64, box.grid_size...)
    slopefn = stencil(Float64, Periodic{2}, (3, 3), slope_kernel)
    slopefn(water_depth, slope, box.phys_scale ./u"m")
    slope_flat = zeros(box.grid_size...)
    slopefn(water_depth_flat, slope_flat, box.phys_scale ./u"m")

    (denudation_mass_HIGH_CO2) = denudation(box, DENUDATION_HIGH_CO2, water_depth, slope,MODEL1,STATE1)
    (denudation_mass_LOW_CO2) = denudation(box, DENUDATION_LOW_CO2, water_depth, slope,MODEL1,STATE1)
    (denudation_mass_LOW_P) = denudation(box, DENUDATION_LOW_P, water_depth, slope,MODEL1,STATE1)
    (denudation_mass_HIGH_P) = denudation(box, DENUDATION_HIGH_P, water_depth, slope,MODEL1,STATE1)

    (denudation_mass_phys) = denudation(box, DENUDATION_PHYS, water_depth, slope,MODEL1,STATE1)
    (denudation_mass_phys_flat) = denudation(box, DENUDATION_PHYS, water_depth_flat, slope_flat,MODEL1,STATE1)
    for idx in CartesianIndices(STATE1.ca)
        f = STATE1.ca[idx]
        if f == 0
            continue
        end

    inf_map[idx] = MODEL1[f].infiltration_coefficient
    end

    # needs restructuring to include the sediment buffer, sediment thickness
    # the resdistribution function should be applied to sediment popped from the buffer
    # not on the denuded mass since this is just the total mass, not split by facies
    # also now that the sediment_thickness is checked, none of these do any denudation!

    @test sum(denudation_mass_HIGH_CO2) > sum(denudation_mass_LOW_CO2)
    @test sum(denudation_mass_LOW_P) < sum(denudation_mass_HIGH_P)
    @test sum(denudation_mass_phys) > sum(denudation_mass_phys_flat)

    (redistribution_mass) = redistribution(box,DENUDATION_PHYS,denudation_mass_phys .*1.0u"Myr",water_depth)
    @test sum(denudation_mass_phys .*1.0u"Myr") ≈ sum(redistribution_mass)

end

@testset "dominant_facies" begin
    state = D.State(
        step=0,
        bathymetry=zeros(Height, 3,3),
        sediment_thickness=zeros(Height, 3,3),
        sediment_buffer = zeros(Float64,10,3,3,3)
    )
    some_sed = zeros(3,3,3)
    some_sed[1,:,:] .+= 1.0 # 1 is extactly enough to tip over into seocnd bucket, so first is empty!
    push_sediment!(state.sediment_buffer, some_sed)
    @test dominant_facies(state, CartesianIndex(1,1)) == 1    

    some_sed[2,:,:] .+= 2.0 
    push_sediment!(state.sediment_buffer, some_sed)
    @test dominant_facies(state, CartesianIndex(1,1)) == 2

    some_sed[3,:,:] .+= 3.0 
    push_sediment!(state.sediment_buffer, some_sed)
    @test dominant_facies(state, CartesianIndex(1,1)) == 3

    # if there's no maximum, takes the first index
    same_sed = ones(3,3,3)
    push_sediment!(state.sediment_buffer, same_sed)
    @test dominant_facies(state, CartesianIndex(1,1)) == 1

end

end