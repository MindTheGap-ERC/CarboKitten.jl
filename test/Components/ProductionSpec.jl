# ~/~ begin <<docs/src/components/production.md#test/Components/ProductionSpec.jl>>[init]
module ProductionSpec
    using Test
    using CarboKitten
    using CarboKitten.Components.Common
    using CarboKitten.Components.Production: Facies, Input, uniform_production
    using CarboKitten.Components.WaterDepth: initial_state
    using CarboKitten.Production: InterpolatedProduction
    # ~/~ begin <<docs/src/components/production.md#production-spec>>[init]
    @testset "Components/Production" begin
        let prod = BenthicProduction(
                maximum_growth_rate = 500u"m/Myr",
                extinction_coefficient = 0.8u"m^-1",
                saturation_intensity = 60u"W/m^2"),
            input = Input(
                box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
                time = TimeProperties(Δt=1.0u"kyr", steps=10),
                sea_level = t -> 0.0u"m",
                initial_topography = (x, y) -> -10u"m",
                subsidence_rate = 0.0u"m/Myr",
                facies = [Facies(production=prod)],
                insolation = 400.0u"W/m^2")
    
            state = initial_state(input)
            prod = uniform_production(input)(state)
            @test all(prod[1:end-1,:] .>= prod[2:end,:])
        end
    
        let prod = PelagicProduction(
                maximum_growth_rate = 5u"1/Myr",
                extinction_coefficient = 0.8u"m^-1",
                saturation_intensity = 60u"W/m^2"),
            input = Input(
                box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
                time = TimeProperties(Δt=1.0u"kyr", steps=10),
                sea_level = t -> 0.0u"m",
                initial_topography = (x, y) -> -10u"m",
                subsidence_rate = 0.0u"m/Myr",
                facies = [Facies(production=prod)],
                insolation = 400.0u"W/m^2")
    
            state = initial_state(input)
            prod = uniform_production(input)(state)
            @test all(prod[1:end-1,:] .<= prod[2:end,:])
        end
    end
    
    @testset "Components/Production/interpolated" begin
        let prod = InterpolatedProduction(
                maximum_production = 500u"m/Myr",
                depth_knots = [0.0u"m", 5.0u"m", 15.0u"m", 50.0u"m"],
                multipliers = [0.0,     1.0,     1.0,      0.0]),
            input = Input(
                box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=5.0u"m"),
                time = TimeProperties(Δt=1.0u"kyr", steps=10),
                sea_level = t -> 0.0u"m",
                initial_topography = (x, y) -> -50u"m",
                subsidence_rate = 0.0u"m/Myr",
                facies = [Facies(production=prod)],
                insolation = 400.0u"W/m^2")
    
            state = initial_state(input)
            p = uniform_production(input)(state)
            # production at 10 m depth (mid-plateau) should be > production at 0 m (multiplier=0)
            @test p[1, 3, 1] > p[1, 1, 1]
            # production should be non-negative everywhere
            @test all(p .>= 0.0u"m")
        end
    end
    # ~/~ end
    # ~/~ begin <<docs/src/components/production.md#production-spec>>[1]
    @testset "Components/Production/variable_insolation" begin
        let prod = BenthicProduction(
                maximum_growth_rate = 500u"m/Myr",
                extinction_coefficient = 0.8u"m^-1",
                saturation_intensity = 60u"W/m^2"),
            input = Input(
                box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
                time = TimeProperties(Δt=1.0u"kyr", steps=10),
                sea_level = t -> 0.0u"m",
                initial_topography = (x, y) -> -10u"m",
                subsidence_rate = 0.0u"m/Myr",
                facies = [Facies(production=prod)],
                insolation = t -> 40.0u"W/m^2/kyr" * t)
    
            state = initial_state(input)
            state.step = 1
            prod1 = copy(uniform_production(input)(state))
            state.step = 10
            prod2 = copy(uniform_production(input)(state))
            @test all(prod2 .> prod1)
        end
    end
    # ~/~ end
end
# ~/~ end
