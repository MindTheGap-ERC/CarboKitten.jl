# ~/~ begin <<docs/src/components/wavefield.md#examples/wavefield/airy_test.jl>>[init]
module AiryTest

using Unitful
using CarboKitten
using CarboKitten.WaveField: AiryWaveField, AiryWaveComponent, energy_flux

# ── 1. Basic sanity: single swell ──────────────────────────────────────

wf = AiryWaveField(components=[
    AiryWaveComponent(amplitude=1.5u"m", period=8.0u"s", direction=0.0),
])

# Evaluate at several depths
for h in [5.0, 10.0, 20.0, 50.0, 100.0]
    v, s = wf(h * u"m")
    E = energy_flux(wf, h * u"m")
    println("h=$(h)m  v=$(v)  shear=$(s)  E=$(E)")
end

# ── 2. Verify breaking: velocity should plateau in very shallow water ──

v_deep, _  = wf(50.0u"m")
v_shallow, _ = wf(0.5u"m")
println("\nDeep velocity: $(v_deep)")
println("Shallow velocity (depth-limited): $(v_shallow)")

# ── 3. Multi-directional sea state ─────────────────────────────────────

wf2 = AiryWaveField(components=[
    AiryWaveComponent(amplitude=1.5u"m", period=8.0u"s", direction=0.0),
    AiryWaveComponent(amplitude=0.5u"m", period=5.0u"s", direction=π/4),
    AiryWaveComponent(amplitude=0.3u"m", period=12.0u"s", direction=-π/6),
])

v, s = wf2(10.0u"m")
println("\nMulti-directional at 10m: v=$(v)  shear=$(s)")

# ── 4. Plug into ALCAP facies (backward compat demo) ───────────────────

input = ALCAP.Input(
    tag = "airy-wave-test",
    box = Box{Coast}(grid_size=(50, 25), phys_scale=150.0u"m"),
    time = TimeProperties(Δt=0.0002u"Myr", steps=1000),
    output = Dict(
        :profile => OutputSpec(slice=(:, 13), write_interval=1)),
    ca_interval = 1,
    initial_topography = (x, y) -> -x / 300.0,
    sea_level = t -> 4.0u"m" * sin(2π * t / 0.2u"Myr"),
    subsidence_rate = 50.0u"m/Myr",
    disintegration_rate = 50.0u"m/Myr",
    lithification_time = 100.0u"yr",
    insolation = 400.0u"W/m^2",
    sediment_buffer_size = 50,
    depositional_resolution = 0.5u"m",
    facies = [ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        production = ALCAP.Example.FACIES[1].production,
        transport_coefficient = 50.0u"m/yr",
        wave_velocity = wf,
        name = "airy-swell")])

mkpath("data/output")
run_model(Model{ALCAP}, input, "data/output/airy-wave-test.h5")
println("\nModel run complete → data/output/airy-wave-test.h5")

end
# ~/~ end
