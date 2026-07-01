module FaciesClassificationSpec

using Test
using Unitful
using GeometryBasics
using CarboKitten
using CarboKitten.FaciesClassification: FaciesRule, classify_block,
    reclassify_data, velocity_magnitude
using CarboKitten.Output.Abstract: Header, Axes, DataSlice, water_depth

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

const TEST_VELOCITY = w -> begin
    v = 1.0u"m/yr" * exp(-w / (20.0u"m"))
    s = -v / (20.0u"m")
    (Vec2(v, 0.0u"m/yr"), Vec2(s, 0.0u"1/yr"))
end

@testset "FaciesClassification/velocity_magnitude sanity" begin
    V_shallow = velocity_magnitude(TEST_VELOCITY, 5.0u"m")
    V_deep    = velocity_magnitude(TEST_VELOCITY, 50.0u"m")
    @test V_shallow > V_deep
    @test unit(V_shallow) == u"m/Myr"
    @test velocity_magnitude(nothing, 5.0u"m") == 0.0u"m/Myr"
end

@testset "FaciesClassification/classify_block — depth gate" begin
    rules = [
        FaciesRule(name="shallow", depth_range=(0.0u"m",  5.0u"m")),
        FaciesRule(name="mid",     depth_range=(5.0u"m",  30.0u"m")),
        FaciesRule(name="deep",    depth_range=(30.0u"m", 200.0u"m")),
    ]
    wv = 0.0u"m/Myr"
    @test classify_block(rules, [1.0], 2.0u"m",   wv) == 1
    @test classify_block(rules, [1.0], 5.0u"m",   wv) == 1
    @test classify_block(rules, [1.0], 15.0u"m",  wv) == 2
    @test classify_block(rules, [1.0], 100.0u"m", wv) == 3
    @test classify_block(rules, [1.0], 500.0u"m", wv) == 4
end

@testset "FaciesClassification/classify_block — fraction gate" begin
    rules = [
        FaciesRule(name="euphotic_dom",
                   sediment_fractions = Dict(1 => (0.6, 1.0)),
                   depth_range = (-Inf*u"m", Inf*u"m")),
        FaciesRule(name="mixed",
                   depth_range = (-Inf*u"m", Inf*u"m")),
    ]
    wv = 0.0u"m/Myr"
    @test classify_block(rules, [0.8, 0.2], 5.0u"m", wv) == 1
    @test classify_block(rules, [0.3, 0.7], 5.0u"m", wv) == 2
    @test classify_block(rules, [0.0, 0.0], 5.0u"m", wv) == 2
end

@testset "FaciesClassification/classify_block — wave velocity gate" begin
    V_high = velocity_magnitude(TEST_VELOCITY, 3.0u"m")
    V_low  = velocity_magnitude(TEST_VELOCITY, 80.0u"m")
    threshold = (V_high + V_low) / 2
    rules = [
        FaciesRule(name="high_velocity", wave_velocity_range=(threshold, Inf*u"m/Myr")),
        FaciesRule(name="low_velocity",  wave_velocity_range=(0.0u"m/Myr", threshold)),
    ]
    @test classify_block(rules, [1.0], 5.0u"m", V_high) == 1
    @test classify_block(rules, [1.0], 5.0u"m", V_low)  == 2
end

@testset "FaciesClassification/classify_block — combined gates" begin
    V_ref = velocity_magnitude(TEST_VELOCITY, 10.0u"m")
    rules = [
        FaciesRule(name="grainstone",
                   sediment_fractions = Dict(1 => (0.5, 1.0)),
                   depth_range          = (0.0u"m", 15.0u"m"),
                   wave_velocity_range = (V_ref, Inf*u"m/Myr")),
        FaciesRule(name="wackestone",
                   depth_range = (0.0u"m", Inf*u"m")),
    ]
    V_hi = velocity_magnitude(TEST_VELOCITY, 5.0u"m")
    V_lo = velocity_magnitude(TEST_VELOCITY, 80.0u"m")
    @test classify_block(rules, [0.7, 0.3], 5.0u"m",  V_hi) == 1
    @test classify_block(rules, [0.7, 0.3], 25.0u"m", V_hi) == 2
    @test classify_block(rules, [0.7, 0.3], 5.0u"m",  V_lo) == 2
    @test classify_block(rules, [0.3, 0.7], 5.0u"m",  V_hi) == 2
end

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
    @test size(new_data.deposition, 1) == 3
    @test size(new_data.deposition, 2) == nx
    @test size(new_data.deposition, 3) == n_t
    @test new_data.sediment_thickness === data.sediment_thickness
    @test new_header.attributes["classified_facies"] == ["A", "B", "fallback"]
end

@testset "FaciesClassification/reclassify_data — bucket routing" begin
    n_f, nx, n_t = 2, 2, 2
    header = _make_header(n_f; nx=nx, n_t=n_t)
    dep = zeros(typeof(1.0u"m"), n_f, nx, n_t)
    dep[1, 1, :] .= 1.0u"m"
    dep[2, 2, :] .= 1.0u"m"
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
    @test all(new_data.deposition[2, 2, :] .≈ 1.0u"m")
end

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
end

@testset "FaciesClassification/reclassify_data — velocity routing" begin
    n_f, nx, n_t = 1, 1, 1
    header = _make_header(n_f; nx=nx, n_t=n_t, sea_level_m=10.0)
    dep  = fill(1.0u"m", n_f, nx, n_t)
    data = _make_slice(n_f, nx, n_t; deposition=dep)
    V_at_10   = velocity_magnitude(TEST_VELOCITY, 10.0u"m")
    threshold = V_at_10 / 2
    rules = [
        FaciesRule(name="wave_active",
                   wave_velocity_range = (threshold, Inf*u"m/Myr"),
                   depth_range = (-Inf*u"m", Inf*u"m")),
        FaciesRule(name="wave_quiet",
                   depth_range = (-Inf*u"m", Inf*u"m")),
    ]
    _, new_data = reclassify_data(header, data, rules; wave_velocity=TEST_VELOCITY)
    @test all(new_data.deposition[1, :, :] .≈ 1.0u"m")
    @test all(new_data.deposition[2, :, :] .≈ 0.0u"m")
end

@testset "FaciesClassification/reclassify_data — no wave velocity fallthrough" begin
    n_f, nx, n_t = 1, 1, 1
    header = _make_header(n_f; nx=nx, n_t=n_t)
    dep  = fill(1.0u"m", n_f, nx, n_t)
    data = _make_slice(n_f, nx, n_t; deposition=dep)
    rules = [
        FaciesRule(name="wave_only",
                   wave_velocity_range = (1.0u"m/Myr", Inf*u"m/Myr"),
                   depth_range = (-Inf*u"m", Inf*u"m")),
        FaciesRule(name="catch_all", depth_range=(-Inf*u"m", Inf*u"m")),
    ]
    _, new_data = reclassify_data(header, data, rules)
    @test all(new_data.deposition[1, :, :] .≈ 0.0u"m")
    @test all(new_data.deposition[2, :, :] .≈ 1.0u"m")
end

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
    @test all(new_data.deposition[1, :, :] .≈ 1.0u"m")
    @test all(new_data.deposition[2, :, :] .≈ 0.0u"m")
end

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

@testset "FaciesClassification/backward-compatibility" begin
    @test !hasfield(CarboKitten.Components.FaciesBase.Facies, :classification_rules)
    @test !hasfield(CarboKitten.Models.ALCAP.Input,           :classification_rules)
    @test !hasfield(CarboKitten.Models.ALCAP.Input,           :save_water_depth)
end

end  # module FaciesClassificationSpec
