module DenudationSpec
using Test
using Unitful

using CarboKitten
using CarboKitten.Denudation.DissolutionMod: dissolution, dominant_facies
using CarboKitten.SedimentStack: push_sediment!, peek_sediment 
using CarboKitten.Components.SedimentBuffer: pop_sediment
using CarboKitten.Components.Common
using CarboKitten.Components: Denudation as D


using CarboKitten.Stencil: Periodic, Reflected, stencil
using CarboKitten.Config: Vectors, TimeProperties
using CarboKitten.Boxes: Box
using CarboKitten.Models: WithDenudation as WD
using CarboKitten.Denudation.EmpiricalDenudationMod: slope_kernel
using CarboKitten.Denudation: denudation, redistribution, Dissolution, NoDenudation, PhysicalErosion, EmpiricalDenudation

FACIES1 = [
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
        erodibility = 0.23u"m/yr",
        initial_sediment=5.0u"m"),

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
        erodibility = 0.23u"m/yr",
        initial_sediment=5.0u"m"),

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
        erodibility = 0.23u"m/yr",
        initial_sediment=5.0u"m"
        )
    ]

function denudation_test_input(denudation_type, sea_level)
    input = WD.Input(
        tag="den_test",
        box=Box{Periodic{2}}(grid_size=(5, 5), phys_scale=1.0u"km"),
        time=TimeProperties(
			Δt = 200.0u"yr",
			steps = 10),
        output=Dict(:profile => OutputSpec(slice=(:, 2), write_interval=1)),
        ca_interval=1,
        initial_topography=(x,y) -> -15.0u"m",
        sea_level=sea_level, # make this into an array and pass as arg?
        facies=FACIES1, # also arg?
        insolation=400.0u"W/m^2",
        denudation=denudation_type
    )
end

@testset "DenudationTST" begin
 
    DENUDATION_HIGH_CO2 = Dissolution(temp = 293.0u"K",precip = 1.0u"m/yr", pco2 = 10^(-1.5)*u"atm",reactionrate = 2e-3u"m/yr")
    DENUDATION_LOW_CO2 = Dissolution(temp = 293.0u"K",precip = 1.0u"m/yr", pco2 = 10^(-2.5)*u"atm",reactionrate = 2e-3u"m/yr")
    DENUDATION_LOW_P = EmpiricalDenudation(precip = 0.8u"m/yr")
    DENUDATION_HIGH_P = EmpiricalDenudation(precip = 1.0u"m/yr")
    DENUDATION_PHYS = PhysicalErosion()
    
    water_depth = -100 .* [ 0.989943  0.48076   0.518983  0.997996   0.895681;
                    0.872733  0.208779  0.882917  0.550494   0.674066;
                    0.57987   0.619433  0.769506  0.593786   0.856186;
                    0.407728  0.469545  0.896348  0.473817   0.797112;
                    0.610194  0.921632  0.322729  0.0103646  0.691191]

    water_depth_flat = -0.5 .* ones(5,5)
    

    # Dissolution
    INPUT_HCO2 = denudation_test_input(DENUDATION_HIGH_CO2, water_depth)

    slope = rand(Float64, 5, 5)
    slopefn = stencil(Float64, Periodic{2}, (3, 3), slope_kernel)
    slopefn(water_depth, slope, INPUT_HCO2.box.phys_scale ./u"m")
    slope_flat = zeros(Float64,5,5)
    slopefn(water_depth_flat, slope_flat, INPUT_HCO2.box.phys_scale ./u"m")

    STATE_HCO2 = WD.initial_state(INPUT_HCO2)
    denudation_mass_HCO2 = denudation(INPUT_HCO2)(STATE_HCO2, water_depth, slope)

    INPUT_LCO2 = denudation_test_input(DENUDATION_LOW_CO2, water_depth)
    STATE_LCO2 = WD.initial_state(INPUT_LCO2)
    denudation_mass_LCO2 = denudation(INPUT_LCO2)(STATE_LCO2, water_depth, slope)

    @test sum(denudation_mass_HCO2) > sum(denudation_mass_LCO2)

    # Empirical 
    INPUT_HP = denudation_test_input(DENUDATION_HIGH_P, water_depth)
    STATE_HP = WD.initial_state(INPUT_HP)
    denudation_mass_HP = denudation(INPUT_HP)(STATE_HP, water_depth, slope)
    
    INPUT_LP = denudation_test_input(DENUDATION_LOW_P, water_depth)
    STATE_LP = WD.initial_state(INPUT_LP)
    denudation_mass_LP = denudation(INPUT_LP)(STATE_LP, water_depth, slope)

    @test sum(denudation_mass_LP) < sum(denudation_mass_HP)

    # physical
    INPUT_PHYS_SLOPE = denudation_test_input(DENUDATION_PHYS, water_depth)
    STATE_PHYS_SLOPE = WD.initial_state(INPUT_PHYS_SLOPE)
    denudation_mass_PHYS_SLOPE = denudation(INPUT_PHYS_SLOPE)(STATE_PHYS_SLOPE, water_depth, slope)
   
    println(denudation_mass_PHYS_SLOPE)

    INPUT_PHYS_FLAT = denudation_test_input(DENUDATION_PHYS, water_depth_flat)
    STATE_PHYS_FLAT = WD.initial_state(INPUT_PHYS_FLAT)
    denudation_mass_PHYS_FLAT = denudation(INPUT_PHYS_FLAT)(STATE_PHYS_FLAT, water_depth_flat, slope_flat)

    @test sum(denudation_mass_PHYS_SLOPE) > sum(denudation_mass_PHYS_FLAT)
    @test sum(denudation_mass_PHYS_FLAT) ≈ 0.0u"m"

    # redistribution
    denuded_sediment = Array{Amount, 3}(undef, 3, INPUT_PHYS_SLOPE.box.grid_size...)
    pop_sediment(INPUT_PHYS_SLOPE)(STATE_PHYS_SLOPE, denudation_mass_PHYS_SLOPE, denuded_sediment)
    redistribution_mass = redistribution(INPUT_PHYS_SLOPE)(STATE_PHYS_SLOPE, water_depth, denuded_sediment)
    
    @test sum(denuded_sediment) ≈ sum(redistribution_mass)

end

@testset "dominant_facies" begin
    state = D.State(
        step=0,
        bathymetry=zeros(Height, 3,3),
        sediment_thickness=zeros(Height, 3,3),
        sediment_buffer = zeros(Float64,10,3,3,3)
    )
    some_sed = zeros(Float64,3,3,3)
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
    same_sed = ones(Float64,3,3,3)
    push_sediment!(state.sediment_buffer, same_sed)
    @test dominant_facies(state, CartesianIndex(1,1)) == 1

end

end