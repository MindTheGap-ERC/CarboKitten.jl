# ~/~ begin <<docs/src/subsidence.md#examples/subsidence.jl>>[init]
# =============================================================================
# Subsidence-modification examples
# =============================================================================
#
# Three runs of the ALCAP example, each varying only the subsidence inputs.
# All three should work side-by-side and produce comparable H5 outputs.
#
# 1. Legacy: uniform scalar subsidence rate (50 m/Myr everywhere).
# 2. Per-cell rate: ramped subsidence increasing along x.
# 3. With modifiers: uniform base rate, but a region is halved over a
#    specific time window and another receives an additive bump.
# =============================================================================

module Script

using Unitful
using CarboKitten
using CarboKitten.Components.WaterDepth: AbstractSubsidenceModifier, MultiplyRate, AddRate, SetRate, Halve

const PATH = "data/output"
const FACIES = ALCAP.Example.FACIES   # reuse the shipped facies definitions

const BOX = Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m")

base_input(tag, subsidence; modifiers=AbstractSubsidenceModifier[]) = ALCAP.Input(
    tag=tag,
    box=BOX,
    time=TimeProperties(Δt=0.0002u"Myr", steps=5000),
    output=Dict(
        :topography => OutputSpec(slice=(:,:), write_interval=10),
        :profile    => OutputSpec(slice=(:, 25), write_interval=1)),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=t -> 4.0u"m" * sin(2π * t / 0.2u"Myr"),
    subsidence_rate=subsidence,
    subsidence_modifiers=modifiers,
    disintegration_rate=50.0u"m/Myr",
    lithification_time=100.0u"yr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    facies=FACIES)

# -- 1. Uniform (legacy path) -------------------------------------------------
function run_scalar()
    input = base_input("subs-scalar", 50.0u"m/Myr")
    run_model(Model{ALCAP}, input, "$(PATH)/subs-scalar.h5")
end

# -- 2. Per-cell rate map -----------------------------------------------------
# Ramp from 30 to 70 m/Myr along x. Build the matrix once and pass it in.
function run_matrix()
    nx, ny = BOX.grid_size
    rates = [30.0u"m/Myr" + 40.0u"m/Myr" * (i - 1) / (nx - 1) for i in 1:nx, _ in 1:ny]
    input = base_input("subs-matrix", rates)
    run_model(Model{ALCAP}, input, "$(PATH)/subs-matrix.h5")
end

# -- 3. Base rate plus localized modifiers -----------------------------------
function run_modifiers()
    modifiers = [
        # Halve subsidence over the first km of x for the first half of the run
        Halve(x_range=(0.0u"m", 1000.0u"m"),
              t_range=(0.0u"Myr", 0.5u"Myr")),
        # Add 20 m/Myr to a central patch, applied for the full run
        AddRate(20.0u"m/Myr";
                x_range=(3.0u"km", 6.0u"km"),
                y_range=(2.0u"km", 5.0u"km")),
        # Pin the rate to 0 over the rightmost strip for the last quarter
        SetRate(0.0u"m/Myr";
                x_range=(12.0u"km", 15.0u"km"),
                t_range=(0.75u"Myr", 1.0u"Myr")),
    ]
    input = base_input("subs-modifiers", 50.0u"m/Myr"; modifiers=modifiers)
    run_model(Model{ALCAP}, input, "$(PATH)/subs-modifiers.h5")
end

function main()
    run_scalar()
    run_matrix()
    run_modifiers()
end

end

Script.main()
# ~/~ end
