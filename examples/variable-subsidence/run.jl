# ~/~ begin <<docs/src/cases/variable-subsidence.md#examples/variable-subsidence/run.jl>>[init]

module VariableSubsidence

using CarboKitten

const PATH = "data/output"
const TAG = "variable-subsidence"

# ~/~ begin <<docs/src/cases/variable-subsidence.md#variable-subsidence-rate>>[init]
const BOX = Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m")
const WIDTH = BOX.grid_size[2] * BOX.phys_scale

const BACKGROUND_RATE = 40.0u"m/Myr"
const FAULT_RATE = 400.0u"m/Myr"
const FAULT_DECAY = 0.25u"Myr"

subsidence(x, y, t) = BACKGROUND_RATE +
    FAULT_RATE * sin(π * y / WIDTH)^2 * exp(-t / FAULT_DECAY)
# ~/~ end

const PROFILES = [:edge => 1, :flank => 13, :centre => 25]

const INPUT = ALCAP.Input(
    tag=TAG,
    box=BOX,
    time=TimeProperties(
        Δt=0.0002u"Myr",
        steps=5000),
    output=Dict(
        :topography => OutputSpec(slice=(:, :), write_interval=100),
        (k => OutputSpec(slice=(:, j), write_interval=1) for (k, j) in PROFILES)...),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=t -> 4.0u"m" * sin(2π * t / 0.2u"Myr"),
    subsidence_rate=SpaceTimeVariableSubsidenceRate(subsidence),
    disintegration_rate=50.0u"m/Myr",
    lithification_time=100.0u"yr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    facies=ALCAP.Example.FACIES)

function main()
    mkpath(PATH)
    run_model(Model{ALCAP}, INPUT, joinpath(PATH, "$(TAG).h5"))
end

end

VariableSubsidence.main()
# ~/~ end
