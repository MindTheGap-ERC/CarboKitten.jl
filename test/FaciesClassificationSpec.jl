# ~/~ begin <<docs/src/facies-classification.md#test/FaciesClassificationSpec.jl>>[init]
module FaciesClassificationSpec

using Test
using Unitful
using CarboKitten
using CarboKitten.FaciesClassification: FaciesRule, classify_block, reclassify_data
using CarboKitten.WaveField: AiryWaveField, AiryWaveComponent, energy_flux
using CarboKitten.Output.Abstract: Header, Axes, Data, DataSlice, water_depth

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function _make_header(n_facies::Int; nx=4, n_t=3,
                      sea_level_m=10.0, subsidence=0.0u"m/Myr")
    t = collect(0.0:0.2:n_t*0.2) * u"Myr"
    Header(
        tag                = "test",
        axes               = Axes(
            x = collect(0.0:150.0:(nx-1)*150.0) * u"m",
            y = [0.0u"m"],
            t = t),
        Δt                 = 0.2u"Myr",
        time_steps         = n_t,
        grid_size          = (nx, 1),
        n_facies           = n_facies,
        initial_topography = zeros(typeof(1.0u"m"), nx, 1),
        sea_level          = fill(sea_level_m * u"m", n_t + 1),
        subsidence_rate    = subsidence,
        data_sets          = Dict(),
        attributes         = Dict())
end

function _make_slice(n_f, nx, n_t;
                     deposition         = zeros(typeof(1.0u"m"), n_f, nx, n_t),
                     sediment_thickness = zeros(typeof(1.0u"m"), nx, n_t),
                     water_depth_arr    = nothing)
    DataSlice(
        slice              = (:, 1),
        write_interval     = 1,
        disintegration     = zeros(typeof(1.0u"m"), n_f, nx, n_t),
        production         = zeros(typeof(1.0u"m"), n_f, nx, n_t),
        deposition         = deposition,
        sediment_thickness = sediment_thickness,
        water_depth        = water_depth_arr)
end

# Minimal single-component wave field used across tests
const TEST_WAVE = AiryWaveField(components=[
    AiryWaveComponent(amplitude=1.5u"m", period=8.0u"s", direction=0.0)])

# ---------------------------------------------------------------------------
# 1.  energy_flux from WaveField — sanity checks
# ---------------------------------------------------------------------------
@testset "FaciesClassification/energy_flux sanity" begin
    # Flux must decrease with depth (energy dissipates toward deep water)
    E_shallow = energy_flux(TEST_WAVE, 5.0u"m")
    E_deep    = energy_flux(TEST_WAVE, 50.0u"m")
    @test E_shallow > E_deep

    # Units must be W/m
    @test unit(E_shallow) == u"W/m"

    # Zero or negative depth: no energy reaches the bed
    @test energy_flux(TEST_WAVE, 0.0u"m") == 0.0u"W/m"
end

# ---------------------------------------------------------------------------
# 2.  classify_block — depth gate
# ---------------------------------------------------------------------------
@testset "FaciesClassification/classify_block — depth gate" begin
    rules = [
        FaciesRule(name="shallow", depth_range=(0.0u"m",  5.0u"m")),
        FaciesRule(name="mid",     depth_range=(5.0u"m",  30.0u"m")),
        FaciesRule(name="deep",    depth_range=(30.0u"m", 200.0u"m")),
    ]
    we = 0.0u"W/m"
    @test classify_block(rules, [1.0], 2.0u"m",   we) == 1   # shallow
    @test classify_block(rules, [1.0], 5.0u"m",   we) == 1   # boundary: inclusive
    @test classify_block(rules, [1.0], 15.0u"m",  we) == 2   # mid
    @test classify_block(rules, [1.0], 100.0u"m", we) == 3   # deep
    @test classify_block(rules, [1.0], 500.0u"m", we) == 4   # fallback
end

# ---------------------------------------------------------------------------
# 3.  classify_block — sediment fraction gate
# ---------------------------------------------------------------------------
@testset "FaciesClassification/classify_block — fraction gate" begin
    rules = [
        FaciesRule(name="euphotic_dom",
                   sediment_fractions = Dict(1 => (0.6, 1.0)),
                   depth_range = (-Inf*u"m", Inf*u"m")),
        FaciesRule(name="mixed",
                   depth_range = (-Inf*u"m", Inf*u"m")),
    ]
    we = 0.0u"W/m"
    @test classify_block(rules, [0.8, 0.2], 5.0u"m", we) == 1
    @test classify_block(rules, [0.3, 0.7], 5.0u"m", we) == 2
    @test classify_block(rules, [0.0, 0.0], 5.0u"m", we) == 2  # empty block
end

# ---------------------------------------------------------------------------
# 4.  classify_block — wave energy gate (W/m)
# ---------------------------------------------------------------------------
@testset "FaciesClassification/classify_block — wave energy gate" begin
    # Get a realistic high-energy and low-energy value from the test wave field
    E_high = energy_flux(TEST_WAVE, 3.0u"m")    # shallow → high flux
    E_low  = energy_flux(TEST_WAVE, 80.0u"m")   # deep    → low flux
    threshold = (E_high + E_low) / 2

    rules = [
        FaciesRule(name="high_energy", wave_energy_range=(threshold, Inf*u"W/m")),
        FaciesRule(name="low_energy",  wave_energy_range=(0.0u"W/m", threshold)),
    ]
    @test classify_block(rules, [1.0], 5.0u"m", E_high) == 1
    @test classify_block(rules, [1.0], 5.0u"m", E_low)  == 2
end

# ---------------------------------------------------------------------------
# 5.  classify_block — combined gates (grainstone scenario)
# ---------------------------------------------------------------------------
@testset "FaciesClassification/classify_block — combined gates" begin
    E_ref = energy_flux(TEST_WAVE, 10.0u"m")
    rules = [
        FaciesRule(name="grainstone",
                   sediment_fractions = Dict(1 => (0.5, 1.0)),
                   depth_range        = (0.0u"m", 15.0u"m"),
                   wave_energy_range  = (E_ref, Inf*u"W/m")),
        FaciesRule(name="wackestone",
                   depth_range = (0.0u"m", Inf*u"m")),
    ]
    E_hi = energy_flux(TEST_WAVE, 5.0u"m")
    E_lo = energy_flux(TEST_WAVE, 80.0u"m")
    @test classify_block(rules, [0.7, 0.3], 5.0u"m",  E_hi) == 1   # all pass
    @test classify_block(rules, [0.7, 0.3], 25.0u"m", E_hi) == 2   # depth fails
    @test classify_block(rules, [0.7, 0.3], 5.0u"m",  E_lo) == 2   # energy fails
    @test classify_block(rules, [0.3, 0.7], 5.0u"m",  E_hi) == 2   # fraction fails
end

# ---------------------------------------------------------------------------
# 6.  reclassify_data — output shape
# ---------------------------------------------------------------------------
@testset "FaciesClassification/reclassify_data — shape" begin
    n_f, nx, n_t = 3, 4, 5
    header = _make_header(n_f; nx=nx, n_t=n_t)
    data   = _make_slice(n_f, nx, n_t)
    rules  = [
        FaciesRule(name="A", depth_range=(-Inf*u"m", Inf*u"m")),
        FaciesRule(name="B", depth_range=(-Inf*u"m", Inf*u"m")),
    ]
    new_header, new_data = reclassify_data(header, data, rules)

    @test new_header.n_facies == 3
    @test size(new_data.deposition,     1) == 3
    @test size(new_data.production,     1) == 3
    @test size(new_data.disintegration, 1) == 3
    @test size(new_data.deposition, 2) == nx
    @test size(new_data.deposition, 3) == n_t
    @test new_data.sediment_thickness === data.sediment_thickness
    @test new_header.attributes["classified_facies"] == ["A", "B", "fallback"]
end

# ---------------------------------------------------------------------------
# 7.  reclassify_data — bucket routing
# ---------------------------------------------------------------------------
@testset "FaciesClassification/reclassify_data — bucket routing" begin
    n_f, nx, n_t = 2, 2, 2
    header = _make_header(n_f; nx=nx, n_t=n_t)
    dep = zeros(typeof(1.0u"m"), n_f, nx, n_t)
    dep[1, 1, :] .= 1.0u"m"   # cell 1: 100% prod-facies 1
    dep[2, 2, :] .= 1.0u"m"   # cell 2: 100% prod-facies 2
    data = _make_slice(n_f, nx, n_t; deposition=dep)

    rules = [
        FaciesRule(name="f1_dom",
                   sediment_fractions = Dict(1 => (0.9, 1.0)),
                   depth_range = (-Inf*u"m", Inf*u"m")),
        FaciesRule(name="f2_dom",
                   sediment_fractions = Dict(2 => (0.9, 1.0)),
                   depth_range = (-Inf*u"m", Inf*u"m")),
    ]
    _, new_data = reclassify_data(header, data, rules)

    @test all(new_data.deposition[1, 1, :] .≈ 1.0u"m")
    @test all(new_data.deposition[2, 1, :] .≈ 0.0u"m")
    @test all(new_data.deposition[3, 1, :] .≈ 0.0u"m")
    @test all(new_data.deposition[1, 2, :] .≈ 0.0u"m")
    @test all(new_data.deposition[2, 2, :] .≈ 1.0u"m")
    @test all(new_data.deposition[3, 2, :] .≈ 0.0u"m")
end

# ---------------------------------------------------------------------------
# 8.  reclassify_data — mass conservation
# ---------------------------------------------------------------------------
@testset "FaciesClassification/reclassify_data — mass conservation" begin
    n_f, nx, n_t = 3, 3, 4
    header = _make_header(n_f; nx=nx, n_t=n_t)
    dep  = rand(typeof(1.0u"m"), n_f, nx, n_t) .* 0.5u"m"
    data = _make_slice(n_f, nx, n_t; deposition=dep)
    rules = [FaciesRule(name="all", depth_range=(-Inf*u"m", Inf*u"m"))]
    _, new_data = reclassify_data(header, data, rules)
    orig_total = dropdims(sum(dep,                 dims=1), dims=1)
    new_total  = dropdims(sum(new_data.deposition, dims=1), dims=1)
    @test all(new_total .≈ orig_total)
    @test all(new_data.deposition[2, :, :] .≈ 0.0u"m")
end

# ---------------------------------------------------------------------------
# 9.  Wave field routing: shallow cell goes to high-energy class
# ---------------------------------------------------------------------------
@testset "FaciesClassification/reclassify_data — AiryWaveField routing" begin
    # Sea level 10 m, zero topo, zero sed → water depth ≈ 10 m everywhere.
    # energy_flux at 10 m is above energy_flux at 50 m.
    n_f, nx, n_t = 1, 1, 1
    header = _make_header(n_f; nx=nx, n_t=n_t, sea_level_m=10.0)
    dep  = fill(1.0u"m", n_f, nx, n_t)
    data = _make_slice(n_f, nx, n_t; deposition=dep)

    E_at_10 = energy_flux(TEST_WAVE, 10.0u"m")
    # Rule fires if energy > half the value at 10 m (guaranteed to match)
    threshold = E_at_10 / 2

    rules = [
        FaciesRule(name="wave_active",
                   wave_energy_range = (threshold, Inf*u"W/m"),
                   depth_range = (-Inf*u"m", Inf*u"m")),
        FaciesRule(name="wave_quiet",
                   depth_range = (-Inf*u"m", Inf*u"m")),
    ]
    _, new_data = reclassify_data(header, data, rules; wave_field=TEST_WAVE)
    @test all(new_data.deposition[1, :, :] .≈ 1.0u"m")   # wave_active
    @test all(new_data.deposition[2, :, :] .≈ 0.0u"m")   # wave_quiet
    @test all(new_data.deposition[3, :, :] .≈ 0.0u"m")   # fallback
end

# ---------------------------------------------------------------------------
# 10.  No wave field: wave_energy_range rules never fire
# ---------------------------------------------------------------------------
@testset "FaciesClassification/reclassify_data — no wave field fallthrough" begin
    n_f, nx, n_t = 1, 1, 1
    header = _make_header(n_f; nx=nx, n_t=n_t)
    dep  = fill(1.0u"m", n_f, nx, n_t)
    data = _make_slice(n_f, nx, n_t; deposition=dep)

    rules = [
        # Requires wave energy > 0; impossible when no wave field supplied
        FaciesRule(name="wave_only",
                   wave_energy_range = (1.0u"W/m", Inf*u"W/m"),
                   depth_range = (-Inf*u"m", Inf*u"m")),
        FaciesRule(name="catch_all", depth_range=(-Inf*u"m", Inf*u"m")),
    ]
    _, new_data = reclassify_data(header, data, rules)   # no wave_field kwarg
    @test all(new_data.deposition[1, :, :] .≈ 0.0u"m")   # wave_only: never fires
    @test all(new_data.deposition[2, :, :] .≈ 1.0u"m")   # catch_all
end

# ---------------------------------------------------------------------------
# 11.  water_depth — stored field takes priority
# ---------------------------------------------------------------------------
@testset "FaciesClassification/water_depth — stored field" begin
    n_f, nx, n_t = 1, 3, 2
    header    = _make_header(n_f; nx=nx, n_t=n_t, sea_level_m=10.0)
    wd_stored = fill(7.0u"m", nx, n_t)
    data_with  = _make_slice(n_f, nx, n_t; water_depth_arr=wd_stored)
    data_plain = _make_slice(n_f, nx, n_t)
    @test water_depth(header, data_with)  === wd_stored
    @test all(water_depth(header, data_plain) .≈ 10.0u"m")
end

@testset "FaciesClassification/reclassify_data — stored depth overrides header" begin
    n_f, nx, n_t = 1, 1, 1
    # Header says sea level = 10 m, but we store depth = 3 m
    header = _make_header(n_f; nx=nx, n_t=n_t, sea_level_m=10.0)
    dep    = fill(1.0u"m", n_f, nx, n_t)
    data   = _make_slice(n_f, nx, n_t;
                         deposition      = dep,
                         water_depth_arr = fill(3.0u"m", nx, n_t))
    rules = [
        FaciesRule(name="shallow", depth_range=(0.0u"m",  5.0u"m")),
        FaciesRule(name="deep",    depth_range=(5.0u"m",  Inf*u"m")),
    ]
    _, new_data = reclassify_data(header, data, rules)
    @test all(new_data.deposition[1, :, :] .≈ 1.0u"m")   # shallow wins
    @test all(new_data.deposition[2, :, :] .≈ 0.0u"m")
end

# ---------------------------------------------------------------------------
# 12.  Fallback class
# ---------------------------------------------------------------------------
@testset "FaciesClassification/reclassify_data — fallback" begin
    n_f, nx, n_t = 1, 2, 1
    header = _make_header(n_f; nx=nx, n_t=n_t)
    dep    = fill(1.0u"m", n_f, nx, n_t)
    data   = _make_slice(n_f, nx, n_t; deposition=dep)
    rules  = [FaciesRule(name="impossible", depth_range=(-100.0u"m", -50.0u"m"))]
    _, new_data = reclassify_data(header, data, rules)
    @test all(new_data.deposition[1, :, :] .≈ 0.0u"m")
    @test all(new_data.deposition[2, :, :] .≈ 1.0u"m")
end

# ---------------------------------------------------------------------------
# 13.  Backward compatibility
# ---------------------------------------------------------------------------
@testset "FaciesClassification/backward-compatibility" begin
    @test !hasfield(CarboKitten.Components.FaciesBase.Facies, :classification_rules)
    @test !hasfield(CarboKitten.Models.ALCAP.Input,           :classification_rules)
    @test !hasfield(CarboKitten.Models.ALCAP.Input,           :save_water_depth)
end

end  # module FaciesClassificationSpec
# ~/~ end
