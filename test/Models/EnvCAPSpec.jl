# ~/~ begin <<docs/src/models/envcap.md#test/Models/EnvCAPSpec.jl>>[init]
module EnvCAPSpec

using Test
using Unitful

using CarboKitten
using CarboKitten.Models: EnvCAP, WithoutCA
using CarboKitten.Models.EnvCAP.EnvMapping: dominant_env, env_to_factory_prior, apply_prior_bias!

const BOX  = Box{Coast}(grid_size=(5, 5), phys_scale=150.0u"m")
const TIME = TimeProperties(Δt=200.0u"yr", steps=10)

const ENV_FACIES = [
    WithoutCA.Facies(
        production = BenthicProduction(
            maximum_growth_rate    = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity   = 60u"W/m^2"),
        transport_coefficient = 10.0u"m/yr"),
    WithoutCA.Facies(
        production = BenthicProduction(
            maximum_growth_rate    = 200u"m/Myr",
            extinction_coefficient = 0.05u"m^-1",
            saturation_intensity   = 60u"W/m^2"),
        transport_coefficient = 5.0u"m/yr"),
]

const STAGE1_INPUT = WithoutCA.Input(
    tag                     = "envcap_stage1_test",
    box                     = BOX,
    time                    = TIME,
    output                  = Dict(:full => OutputSpec(write_interval=1)),
    initial_topography      = (x, y) -> -x / 300.0,
    insolation              = 400.0u"W/m^2",
    subsidence_rate         = 50.0u"m/Myr",
    disintegration_rate     = 20.0u"m/Myr",
    lithification_time      = 100.0u"yr",
    sediment_buffer_size    = 10,
    depositional_resolution = 0.5u"m",
    facies = ENV_FACIES)

const FACTORY_FACIES = [
    EnvCAP.Facies(
        viability_range  = (4, 10),
        activation_range = (6, 10),
        production = BenthicProduction(
            maximum_growth_rate    = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity   = 60u"W/m^2"),
        transport_coefficient = 10.0u"m/yr"),
    EnvCAP.Facies(
        viability_range  = (4, 10),
        activation_range = (6, 10),
        production = BenthicProduction(
            maximum_growth_rate    = 200u"m/Myr",
            extinction_coefficient = 0.05u"m^-1",
            saturation_intensity   = 60u"W/m^2"),
        transport_coefficient = 5.0u"m/yr"),
]

@testset "EnvMapping/dominant_env" begin
    dep = zeros(typeof(1.0u"m"), 2, 2, 2, 3)
    dep[1, 1, 1, :] .= 1.0u"m"
    dep[2, 2, 2, :] .= 2.0u"m"
    ef = dominant_env(dep)
    @test size(ef) == (2, 2)
    @test ef[1, 1] == 1
    @test ef[2, 2] == 2
    @test ef[1, 2] == 0
end

@testset "EnvMapping/env_to_factory_prior" begin
    env_field = [1 2; 0 1]
    mapping = [0.8 0.2; 0.3 0.7]
    prior = env_to_factory_prior(env_field, mapping)
    @test size(prior) == (2, 2, 2)
    @test prior[1, 1, 1] ≈ 0.8
    @test prior[2, 2, 1] ≈ 0.7
    @test all(prior[:, 1, 2] .≈ 0.5)
end

@testset "EnvMapping/apply_prior_bias! α=0" begin
    import Random: MersenneTwister
    ca    = [1 2; 0 1]
    prior = ones(Float64, 2, 2, 2) .* 0.5
    ca_before = copy(ca)
    apply_prior_bias!(ca, prior, 0.0, MersenneTwister(0))
    @test ca == ca_before
end

@testset "EnvMapping/apply_prior_bias! α=1 zero-prob" begin
    import Random: MersenneTwister
    ca    = ones(Int, 3, 3)
    prior = zeros(Float64, 2, 3, 3)
    prior[2, :, :] .= 1.0
    apply_prior_bias!(ca, prior, 1.0, MersenneTwister(0))
    @test all(ca .== 0)
end

const STAGE1_OUT = run_model(Model{WithoutCA}, STAGE1_INPUT, MemoryOutput(STAGE1_INPUT))

@testset "Stage-1 WithoutCA output shape" begin
    dep = STAGE1_OUT.data_volumes[:full].deposition
    @test size(dep, 1) == 2
    @test size(dep, 2) == 5
    @test size(dep, 3) == 5
end

const ENV_FIELD     = dominant_env(STAGE1_OUT.data_volumes[:full].deposition)
const MAPPING       = [0.9 0.1; 0.1 0.9]
const FACTORY_PRIOR = env_to_factory_prior(ENV_FIELD, MAPPING)

@testset "Factory prior shape and validity" begin
    @test size(FACTORY_PRIOR) == (2, 5, 5)
    @test all(0.0 .<= FACTORY_PRIOR .<= 1.0)
    col_sums = dropdims(sum(FACTORY_PRIOR; dims=1), dims=1)
    @test all(isapprox.(col_sums, 1.0; atol=1e-10))
end

function make_input(ca_refinement)
    EnvCAP.Input(
        tag                     = "envcap_test",
        box                     = BOX,
        time                    = TIME,
        output                  = Dict(:full => OutputSpec(write_interval=1)),
        initial_topography      = (x, y) -> -x / 300.0,
        insolation              = 400.0u"W/m^2",
        subsidence_rate         = 50.0u"m/Myr",
        disintegration_rate     = 20.0u"m/Myr",
        lithification_time      = 100.0u"yr",
        sediment_buffer_size    = 10,
        depositional_resolution = 0.5u"m",
        ca_interval             = 1,
        ca_refinement           = ca_refinement,
        factory_prior           = FACTORY_PRIOR,
        facies = FACTORY_FACIES)
end

const OUT_0  = run_model(Model{EnvCAP}, make_input(0.0), MemoryOutput(make_input(0.0)))
const OUT_05 = run_model(Model{EnvCAP}, make_input(0.5), MemoryOutput(make_input(0.5)))

@testset "Models/EnvCAP output shape" begin
    @test all(size(OUT_0.data_volumes[:full].sediment_thickness)  .== (5, 5, 11))
    @test all(size(OUT_05.data_volumes[:full].sediment_thickness) .== (5, 5, 11))
end

@testset "Models/EnvCAP ca_refinement=0 produces deposition" begin
    @test any(OUT_0.data_volumes[:full].deposition .> 0.0u"m")
end

@testset "Models/EnvCAP factory_prior=nothing behaves like ALCAP" begin
    input_no_prior = EnvCAP.Input(
        tag                     = "envcap_no_prior",
        box                     = BOX,
        time                    = TIME,
        output                  = Dict(:full => OutputSpec(write_interval=1)),
        initial_topography      = (x, y) -> -x / 300.0,
        insolation              = 400.0u"W/m^2",
        subsidence_rate         = 50.0u"m/Myr",
        disintegration_rate     = 20.0u"m/Myr",
        lithification_time      = 100.0u"yr",
        sediment_buffer_size    = 10,
        depositional_resolution = 0.5u"m",
        ca_interval             = 1,
        factory_prior           = nothing,
        ca_refinement           = 1.0,
        facies = FACTORY_FACIES)
    out = run_model(Model{EnvCAP}, input_no_prior, MemoryOutput(input_no_prior))
    @test any(out.data_volumes[:full].deposition .> 0.0u"m")
end

end
# ~/~ end
