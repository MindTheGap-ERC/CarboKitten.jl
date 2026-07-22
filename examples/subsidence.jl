# ~/~ begin <<docs/src/components/waterdepth.md#examples/subsidence.jl>>[init]
module Script

using Unitful
using CarboKitten
using CarboKitten.Components.WaterDepth: MultiplyRate, AddRate, SetRate, Halve,
    AbstractSubsidenceModifier

const PATH   = "data/output"
const FACIES = ALCAP.Example.FACIES
const BOX    = Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m")

function base_input(tag, subsidence_rate; modifiers=AbstractSubsidenceModifier[])
    ALCAP.Input(
        tag  = tag,
        box  = BOX,
        time = TimeProperties(Δt=0.0002u"Myr", steps=5000),
        output = Dict(
            :topography => OutputSpec(slice=(:,:), write_interval=10),
            :profile    => OutputSpec(slice=(:, 25), write_interval=1)),
        ca_interval          = 1,
        initial_topography   = (x, y) -> -x / 300.0,
        sea_level            = t -> 4.0u"m" * sin(2π * t / 0.2u"Myr"),
        subsidence_rate      = subsidence_rate,
        subsidence_modifiers = modifiers,
        disintegration_rate  = 50.0u"m/Myr",
        lithification_time   = 100.0u"yr",
        insolation           = 400.0u"W/m^2",
        sediment_buffer_size = 50,
        depositional_resolution = 0.5u"m",
        facies = FACIES)
end

# 1. Scalar — legacy path, bit-identical to main branch
function run_scalar()
    run_model(Model{ALCAP}, base_input("subs-scalar", 50.0u"m/Myr"),
              "$(PATH)/subs-scalar.h5")
end

# 2. Per-cell rate map — ramp from 30 to 70 m/Myr along x
function run_matrix()
    nx, ny = BOX.grid_size
    rates  = [30.0u"m/Myr" + 40.0u"m/Myr" * (i - 1) / (nx - 1)
              for i in 1:nx, _ in 1:ny]
    run_model(Model{ALCAP}, base_input("subs-matrix", rates),
              "$(PATH)/subs-matrix.h5")
end

# 3. Uniform base rate + localized modifiers
function run_modifiers()
    mods = [
        Halve(x_range=(0.0u"m",     1500.0u"m"),
              t_range=(0.0u"Myr",   0.5u"Myr")),
        AddRate(20.0u"m/Myr";
                x_range=(4500.0u"m",  7500.0u"m"),
                y_range=(3000.0u"m",  6000.0u"m")),
        SetRate(0.0u"m/Myr";
                x_range=(13500.0u"m", 15000.0u"m"),
                t_range=(0.75u"Myr",  1.0u"Myr")),
    ]
    run_model(Model{ALCAP}, base_input("subs-modifiers", 50.0u"m/Myr"; modifiers=mods),
              "$(PATH)/subs-modifiers.h5")
end

function main()
    mkpath(PATH)
    run_scalar()
    run_matrix()
    run_modifiers()
end

end

Script.main()
# ~/~ end
