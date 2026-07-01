# Post-Deposition Facies Classification

CarboKitten models record deposition in terms of **production facies**. This
module reclassifies those production facies into depositional facies during
post-processing.

This updated version removes the dependency on the experimental `WaveField`
module. Wave exposure is now evaluated with the standard CarboKitten
`wave_velocity` formulation:

```julia
wave_velocity(depth) -> (velocity, shear)
```

Only the magnitude of `velocity` is used by the facies rules.

## Design overview

```text
run_model
  -> Data production facies
  -> reclassify_data(rules; wave_velocity=...)
  -> Data classified facies
```

## Main API change

Old Airy-wave version:

```julia
FaciesRule(
    name = "grainstone",
    wave_energy_range = (500.0u"W/m", Inf*u"W/m"))

new_header, new_vol = reclassify_data(header, vol, rules; wave_field=wf)
```

New velocity-based version:

```julia
FaciesRule(
    name = "grainstone",
    wave_velocity_range = (0.20u"m/yr", Inf*u"m/yr"))

new_header, new_vol = reclassify_data(header, vol, rules; wave_velocity=wave_velocity)
```

## User-facing example

```julia
using CarboKitten
using CarboKitten.Export: read_volume
using CarboKitten.FaciesClassification: FaciesRule, reclassify_data
using GeometryBasics
using Unitful

v_prof(v_max, max_depth) = w -> begin
    k = sqrt(0.5) / max_depth
    A = 3.331 * v_max
    α = tanh(k * w)
    β = exp(-k * w)
    v = A * α * β
    s = -A * k * β * (1 - α - α^2)
    (Vec2(v, 0.0u"m/yr"), Vec2(s, 0.0u"1/yr"))
end

wave_velocity = v_prof(0.5u"m/yr", 20.0u"m")

header, vol = read_volume("output/run.h5", :topography)

rules = [
    FaciesRule(
        name                = "grainstone",
        sediment_fractions  = Dict(1 => (0.5, 1.0)),
        depth_range         = (0.0u"m", 20.0u"m"),
        wave_velocity_range = (0.20u"m/yr", Inf*u"m/yr")),
    FaciesRule(
        name                = "packstone",
        sediment_fractions  = Dict(1 => (0.2, 0.8)),
        depth_range         = (0.0u"m", 40.0u"m")),
    FaciesRule(
        name                = "mudstone",
        depth_range         = (10.0u"m", Inf*u"m")),
]

new_header, new_vol = reclassify_data(header, vol, rules; wave_velocity=wave_velocity)
```

## Source file

```{.julia file=src/FaciesClassification.jl}
"""
    FaciesClassification

Post-processing tools for reclassifying CarboKitten production facies into
depositional facies.

Classification is based on:
- production-facies fractions,
- palaeo-water depth,
- optional wave-velocity magnitude.

The wave-velocity formulation follows the standard CarboKitten transport
interface: a `wave_velocity` function takes water depth and returns
`(velocity, shear)`, where velocity has units of length/time and shear has
units of 1/time. Classification uses only the velocity magnitude.
"""
module FaciesClassification

using Unitful
using HDF5
using ..Output.Abstract: Data, Header, water_depth

export FaciesRule, classify_block, reclassify_data, reclassify_volume,
       save_classified, velocity_magnitude

"""
    FaciesRule(; name, sediment_fractions, depth_range, wave_velocity_range)

A single rule for post-depositional facies classification.

## Fields
- `name::String` — label assigned to blocks that match this rule.
- `sediment_fractions::Union{Nothing, Dict{Int,NTuple{2,Float64}}}` —
  optional per-production-facies fraction constraints. Keys are 1-based
  production-facies indices; values are `(lo, hi)` with `0 ≤ lo ≤ hi ≤ 1`.
  All listed constraints must be satisfied simultaneously. `nothing` means
  no constraint on sediment proportion.
- `depth_range::NTuple{2,<:Quantity}` — `(min_depth, max_depth)` water-depth
  window in metres. Bounds are inclusive. Use `(-Inf*u"m", Inf*u"m")` to
  leave depth unconstrained.
- `wave_velocity_range::NTuple{2,<:Quantity}` — `(lo, hi)` for the magnitude
  of the CarboKitten `wave_velocity(depth)[1]` vector. Use
  `(0.0u"m/Myr", Inf*u"m/Myr")` to leave wave velocity unconstrained.

Rules are evaluated in order; first match wins. Unmatched blocks go to the
fallback class at index `length(rules) + 1`.
"""
Base.@kwdef struct FaciesRule
    name::String
    sediment_fractions::Union{Nothing,Dict{Int,NTuple{2,Float64}}} = nothing
    depth_range::NTuple{2,<:Quantity}         = (-Inf * u"m",     Inf * u"m")
    wave_velocity_range::NTuple{2,<:Quantity} = (0.0  * u"m/Myr", Inf * u"m/Myr")
end

"""
    velocity_magnitude(wave_velocity, depth) -> Quantity

Evaluate a CarboKitten-style wave-velocity function at a water depth and return
the magnitude of the velocity vector.

The expected interface is the standard active-layer transport interface:

```julia
wave_velocity(depth) -> (velocity, shear)
```

where `velocity` is a scalar, vector, tuple, or `Vec2` with velocity units.
When `wave_velocity === nothing`, the returned value is `0 m/Myr`.
"""
function velocity_magnitude(wave_velocity, depth::Quantity)
    wave_velocity === nothing && return 0.0u"m/Myr"

    v, _ = wave_velocity(depth)

    if v isa Quantity
        return abs(uconvert(u"m/Myr", v))
    end

    mag = sqrt(sum(abs2, v))
    return uconvert(u"m/Myr", mag)
end

"""
    classify_block(rules, fractions, depth, wave_velocity) -> Int

Classify a single deposited block against an ordered list of `FaciesRule`s.

Returns the 1-based index of the first matching rule, or `length(rules) + 1`
(fallback) when no rule matches.
"""
function classify_block(rules::AbstractVector{FaciesRule},
                        fractions::AbstractVector{<:Real},
                        depth::Quantity,
                        wave_velocity::Quantity)::Int
    for (i, rule) in enumerate(rules)
        depth < rule.depth_range[1] && continue
        depth > rule.depth_range[2] && continue
        wave_velocity < rule.wave_velocity_range[1] && continue
        wave_velocity > rule.wave_velocity_range[2] && continue

        if rule.sediment_fractions !== nothing
            ok = true
            for (fid, (lo, hi)) in rule.sediment_fractions
                frac = length(fractions) >= fid ? fractions[fid] : 0.0
                if !(lo <= frac <= hi)
                    ok = false
                    break
                end
            end
            ok || continue
        end

        return i
    end

    return length(rules) + 1
end

"""
    reclassify_data(header, data, rules; wave_velocity=nothing) -> (Header, Data)

Reclassify a CarboKitten `Data` object using an ordered vector of `FaciesRule`s.

`wave_velocity` is optional. If provided, it must follow the standard
CarboKitten active-layer transport interface:

```julia
wave_velocity(depth) -> (velocity, shear)
```

Only the magnitude of `velocity` is used for classification. This keeps
classification consistent with the original CarboKitten velocity-based
transport formulation and avoids depending on an external `WaveField` module.

Water depth is read from `data.water_depth` when stored, otherwise reconstructed
from the header.

Returns `(new_header, new_data)` with the same `Data{F,D}` type, fully
compatible with existing visualisation routines.
"""
function reclassify_data(header::Header,
                         data::Data{F,D},
                         rules::AbstractVector{FaciesRule};
                         wave_velocity=nothing) where {F,D}
    n_prod  = size(data.deposition, 1)
    n_class = length(rules) + 1
    dep_sz  = size(data.deposition)
    sp_size = dep_sz[2:end-1]
    n_t     = dep_sz[end]

    zero_like(a) = zeros(eltype(a), n_class, sp_size..., n_t)
    dep_out  = zero_like(data.deposition)
    prod_out = zero_like(data.production)
    dis_out  = zero_like(data.disintegration)

    wd_array = water_depth(header, data)

    _wd(sp_idx::CartesianIndex{0}, t) = wd_array[t]
    _wd(sp_idx, t)                    = wd_array[sp_idx, t]

    for t_idx in 1:n_t
        for sp_idx in CartesianIndices(sp_size)
            wd_val = _wd(sp_idx, t_idx)

            dep_col  = @view data.deposition[:,     sp_idx, t_idx]
            prod_col = @view data.production[:,     sp_idx, t_idx]
            dis_col  = @view data.disintegration[:, sp_idx, t_idx]

            total_dep = sum(dep_col)
            fractions = if total_dep > zero(eltype(dep_col))
                ustrip.(dep_col ./ total_dep)
            else
                zeros(Float64, n_prod)
            end

            wv  = velocity_magnitude(wave_velocity, wd_val)
            cls = classify_block(rules, fractions, wd_val, wv)

            dep_out[cls,  sp_idx, t_idx] += total_dep
            prod_out[cls, sp_idx, t_idx] += sum(prod_col)
            dis_out[cls,  sp_idx, t_idx] += sum(dis_col)
        end
    end

    new_header = Header(
        tag                = header.tag,
        axes               = header.axes,
        Δt                 = header.Δt,
        time_steps         = header.time_steps,
        grid_size          = header.grid_size,
        n_facies           = n_class,
        initial_topography = header.initial_topography,
        sea_level          = header.sea_level,
        subsidence_rate    = header.subsidence_rate,
        data_sets          = header.data_sets,
        attributes         = merge(header.attributes,
                                   Dict("facies_classification" => true,
                                        "classified_facies" =>
                                            [[r.name for r in rules]; "fallback"])))

    new_data = Data{F,D}(
        slice              = data.slice,
        write_interval     = data.write_interval,
        disintegration     = dis_out,
        production         = prod_out,
        deposition         = dep_out,
        sediment_thickness = data.sediment_thickness,
        active_layer       = nothing)

    return new_header, new_data
end

"""
    reclassify_volume(header, vol, rules; wave_velocity=nothing) -> (Header, DataVolume)

Convenience wrapper: reclassify a full `DataVolume` into depositional
environments. Identical to calling `reclassify_data` on the volume directly.
"""
function reclassify_volume(header::Header,
                           vol,
                           rules::AbstractVector{FaciesRule};
                           wave_velocity=nothing)
    return reclassify_data(header, vol, rules; wave_velocity=wave_velocity)
end

"""
    save_classified(filename, header, data; group=:classified)

Write a classified `Data` object to HDF5 in the same layout as a standard
CarboKitten output file, so it can be read back with `read_volume`,
`read_slice`, or `read_column`.
"""
function save_classified(filename::AbstractString,
                         header::Header,
                         data::Data;
                         group::Symbol = :classified)

    dep_f    = ustrip.(u"m", data.deposition)
    prod_f   = ustrip.(u"m", data.production)
    disint_f = ustrip.(u"m", data.disintegration)
    thick_f  = ustrip.(u"m", data.sediment_thickness)

    n_facies = header.n_facies
    x_m      = ustrip.(u"m", header.axes.x)
    y_m      = ustrip.(u"m", header.axes.y)
    t_myr    = ustrip.(u"Myr", header.axes.t)
    sl_m     = ustrip.(u"m", header.sea_level)
    topo_m   = ustrip.(u"m", header.initial_topography)

    slice_str(::Colon) = ":"
    slice_str(i::Int)  = string(i)

    h5open(filename, "w") do fid
        grp_in = create_group(fid, "input")
        grp_in["x"]                  = x_m
        grp_in["y"]                  = y_m
        grp_in["t"]                  = t_myr
        grp_in["initial_topography"] = topo_m
        grp_in["sea_level"]          = sl_m

        a = HDF5.attributes(grp_in)
        a["tag"]             = header.tag
        a["delta_t"]         = ustrip(u"Myr", header.Δt)
        a["time_steps"]      = header.time_steps
        a["n_facies"]        = n_facies
        a["subsidence_rate"] = ustrip(u"m/Myr", header.subsidence_rate)

        if haskey(header.attributes, "classified_facies")
            a["classified_facies"] = join(header.attributes["classified_facies"], ",")
        end

        grp = create_group(fid, string(group))
        ag  = HDF5.attributes(grp)
        ag["write_interval"] = data.write_interval
        ag["slice"]          = join(slice_str.(data.slice), ",")

        grp["deposition"]         = dep_f
        grp["production"]         = prod_f
        grp["disintegration"]     = disint_f
        grp["sediment_thickness"] = thick_f
    end

    @info "Classified output saved -> $filename  (group=$group, n_facies=$n_facies)"
    return filename
end

end  # module FaciesClassification

```

## Test specification

```{.julia file=test/FaciesClassificationSpec.jl}
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

```
