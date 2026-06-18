# ~/~ begin <<docs/src/models/envcap.md#examples/model/envcap/run.jl>>[init]

module Script
using CarboKitten
using CarboKitten.Models: EnvCAP, WithoutCA
using CarboKitten.Models.EnvCAP.EnvMapping: dominant_env_block, env_to_factory_prior_block
using CarboKitten.Export: read_volume
using HDF5

# ── shared grid / time ────────────────────────────────────────────────────────
const BOX  = Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m")
const STEPS = 5000
const ΔT    = 200.0u"yr"

# ── Stage 1: WithoutCA environments ──────────────────────────────────────────
#
# Three large-scale environments defined by production curves:
#   Env 1 (inner platform) – high production, steep extinction
#   Env 2 (shoal/margin)   – medium production, moderate extinction
#   Env 3 (slope)          – low production, gentle extinction
const ENV_FACIES = [
    WithoutCA.Facies(
        name = "inner_platform",
        production = BenthicProduction(
            maximum_growth_rate    = 600u"m/Myr",
            extinction_coefficient = 0.9u"m^-1",
            saturation_intensity   = 60u"W/m^2"),
        transport_coefficient = 15.0u"m/yr"),
    WithoutCA.Facies(
        name = "shoal",
        production = BenthicProduction(
            maximum_growth_rate    = 350u"m/Myr",
            extinction_coefficient = 0.3u"m^-1",
            saturation_intensity   = 60u"W/m^2"),
        transport_coefficient = 8.0u"m/yr"),
    WithoutCA.Facies(
        name = "slope",
        production = BenthicProduction(
            maximum_growth_rate    = 120u"m/Myr",
            extinction_coefficient = 0.02u"m^-1",
            saturation_intensity   = 60u"W/m^2"),
        transport_coefficient = 4.0u"m/yr"),
]

const STAGE1_INPUT = WithoutCA.Input(
    tag                     = "envcap_stage1",
    box                     = BOX,
    time                    = TimeProperties(Δt=ΔT, steps=STEPS),
    output                  = Dict(:full => OutputSpec(write_interval=50)),
    initial_topography      = (x, y) -> -x / 300.0,
    sea_level               = t -> 4.0u"m" * sin(2π * t / 200.0u"kyr"),
    subsidence_rate         = 50.0u"m/Myr",
    insolation              = 400.0u"W/m^2",
    disintegration_rate     = 20.0u"m/Myr",
    lithification_time      = 100.0u"yr",
    sediment_buffer_size    = 50,
    depositional_resolution = 0.5u"m",
    facies = ENV_FACIES)

# ── Stage 2: EnvCAP factory facies ───────────────────────────────────────────
#
# Three carbonate factories:
#   Factory 1 (euphotic)   – shallow-water high-energy
#   Factory 2 (oligophotic)– intermediate depth
#   Factory 3 (aphotic)    – deeper, lower energy
const FACTORY_FACIES = [
    EnvCAP.Facies(
        name             = "euphotic",
        viability_range  = (4, 10),
        activation_range = (6, 10),
        production = BenthicProduction(
            maximum_growth_rate    = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity   = 60u"W/m^2"),
        transport_coefficient = 50.0u"m/yr"),
    EnvCAP.Facies(
        name             = "oligophotic",
        viability_range  = (4, 10),
        activation_range = (6, 10),
        production = BenthicProduction(
            maximum_growth_rate    = 400u"m/Myr",
            extinction_coefficient = 0.1u"m^-1",
            saturation_intensity   = 60u"W/m^2"),
        transport_coefficient = 25.0u"m/yr"),
    EnvCAP.Facies(
        name             = "aphotic",
        viability_range  = (4, 10),
        activation_range = (6, 10),
        production = BenthicProduction(
            maximum_growth_rate    = 100u"m/Myr",
            extinction_coefficient = 0.005u"m^-1",
            saturation_intensity   = 60u"W/m^2"),
        transport_coefficient = 12.5u"m/yr"),
]

# Environment-to-factory mapping (3 envs × 3 factories):
#   Inner platform → mostly euphotic
#   Shoal          → euphotic / oligophotic mix
#   Slope          → mostly aphotic
const ENV_TO_FACTORY = [
    0.8  0.15  0.05;   # inner platform
    0.3  0.55  0.15;   # shoal
    0.05 0.25  0.70    # slope
]

function make_stage2_input(factory_prior, ca_refinement; tag)
    EnvCAP.Input(
        tag                     = tag,
        box                     = BOX,
        time                    = TimeProperties(Δt=ΔT, steps=STEPS),
        output                  = Dict(
            :topography => OutputSpec(write_interval=50),
            :profile    => OutputSpec(slice=(:, 25)),
            :full       => OutputSpec(write_interval=50)),
        initial_topography      = (x, y) -> -x / 300.0,
        sea_level               = t -> 4.0u"m" * sin(2π * t / 200.0u"kyr"),
        subsidence_rate         = 50.0u"m/Myr",
        insolation              = 400.0u"W/m^2",
        ca_interval             = 1,
        disintegration_rate     = 20.0u"m/Myr",
        lithification_time      = 100.0u"yr",
        sediment_buffer_size    = 50,
        depositional_resolution = 0.5u"m",
        factory_prior           = factory_prior,
        ca_refinement           = ca_refinement,
        env_random_seed         = 42,
        facies = FACTORY_FACIES)
end

function main()
    # ── Stage 1 ──
    run_model(Model{WithoutCA}, STAGE1_INPUT, "data/output/envcap_stage1.h5")

    # Extract environmental belt from stage-1 output
    header, stage1_data = read_volume("data/output/envcap_stage1.h5", :full)
        env_belt = dominant_env_block(
            stage1_data.deposition;
            dz = 0.5,
            min_nz = 41,
        )
        
        factory_prior = env_to_factory_prior_block(
            env_belt,
            ENV_TO_FACTORY,
        )

    # Save the env_field and factory_prior as attributes in each stage-2 file
    # (they are inputs, not outputs, but useful to archive alongside results)

    # ── Stage 2: three refinement levels, same prior ──
    for (α, tag, path) in [
            (0.0, "envcap_ref0",  "data/output/envcap_ref0.h5"),
            (0.5, "envcap_ref05", "data/output/envcap_ref05.h5"),
            (1.0, "envcap_ref1",  "data/output/envcap_ref1.h5")]
        input = make_stage2_input(factory_prior, α; tag=tag)
        run_model(Model{EnvCAP}, input, path)
    end
end
end

Script.main()
# ~/~ end
