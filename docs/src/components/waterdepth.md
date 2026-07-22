# Water Depth

```component-dag
CarboKitten.Components.WaterDepth
```

The `WaterDepth` module computes the water depth, given the bedrock elevation, sea level curve, subsidence rate and current sediment height.

## Input

- `initial_topography(x, y)` (a.k.a. initial depth) should be a function taking two coordinates in units of meters, returning an elevation also in meters.
- `sea_level(t)` should be a function taking a time in millions of years (Myr) returning the eustatic sealevel. This could also be an interpolated table.
- `subsidence_rate` a constant rate of subsidence in m/Myr. Optionally a `Matrix{Rate}` or `(x, y) -> Rate` function for spatially varying subsidence.
- `subsidence_modifiers` an optional list of modifiers that locally alter the rate in (x, y, t) boxes.

The signs of these quantities should be such that the following equation holds:

$$T + E = S + W,$$

saying Tectonic subsidence plus Eustatic sea-level change equals Sedimentation plus change in Water depth.

## Spatio-temporal subsidence

`subsidence_rate` accepts a scalar `Quantity` (legacy, fast path), a `Matrix{Rate}`, or a function `(x, y) -> Rate`. The optional `subsidence_modifiers` list allows piecewise modification of the rate inside (x, y, t) boxes. All existing scalar inputs with no modifiers are bit-identical to the previous behaviour.

``` {.julia #waterdepth-input}
@kwdef struct Input <: AbstractInput
    sea_level = t -> 0.0u"m"
    initial_topography = (x, y) -> 0.0u"m"
    # scalar Rate (legacy), Matrix{Rate}, or (x,y)->Rate function
    subsidence_rate = 0.0u"m/Myr"
    subsidence_modifiers::Vector{AbstractSubsidenceModifier} = AbstractSubsidenceModifier[]
end
```

`subsidence_rate_map` normalises any of the three forms to a `Matrix{Rate}`:

``` {.julia #subsidence-rate-map}
"""
    subsidence_rate_map(input) -> Matrix{Rate}

Normalise `input.subsidence_rate` to a `Matrix{Rate}` of size
`input.box.grid_size`. Accepts a scalar `Quantity`, a `Matrix{Rate}`, or a
function `(x, y) -> Rate`.
"""
function subsidence_rate_map(input::AbstractInput)
    sr = input.subsidence_rate
    if sr isa AbstractMatrix
        @assert size(sr) == input.box.grid_size
        return sr
    elseif sr isa Function
        x, y = box_axes(input.box)
        return sr.(x, y')
    else
        return fill(sr, input.box.grid_size...)
    end
end
```

## Modifier types

A modifier alters the subsidence rate inside an (x, y, t) box. Ranges accept
`:` (no restriction) or a `(lo, hi)` tuple of `Quantity`.

``` {.julia #subsidence-modifiers}
@kwdef struct MultiplyRate <: AbstractSubsidenceModifier
    factor::Float64
    x_range::AxisSpec = (:)
    y_range::AxisSpec = (:)
    t_range::TimeSpec = (:)
end
MultiplyRate(factor::Real; kwargs...) = MultiplyRate(; factor=Float64(factor), kwargs...)

Halve(;  kwargs...) = MultiplyRate(0.5; kwargs...)
Double(; kwargs...) = MultiplyRate(2.0; kwargs...)

@kwdef struct AddRate <: AbstractSubsidenceModifier
    delta::Rate
    x_range::AxisSpec = (:)
    y_range::AxisSpec = (:)
    t_range::TimeSpec = (:)
end
AddRate(delta::Rate; kwargs...) = AddRate(; delta=delta, kwargs...)

@kwdef struct SetRate <: AbstractSubsidenceModifier
    rate::Rate
    x_range::AxisSpec = (:)
    y_range::AxisSpec = (:)
    t_range::TimeSpec = (:)
end
SetRate(rate::Rate; kwargs...) = SetRate(; rate=rate, kwargs...)

apply_rate(m::MultiplyRate, r::Rate) = m.factor * r
apply_rate(m::AddRate,      r::Rate) = r + m.delta
apply_rate(m::SetRate,      _::Rate) = m.rate
```

| Type           | Effect                         | Key field         |
|----------------|--------------------------------|-------------------|
| `MultiplyRate` | Multiply rate by a factor      | `factor::Float64` |
| `AddRate`      | Add a constant to the rate     | `delta::Rate`     |
| `SetRate`      | Override rate to a fixed value | `rate::Rate`      |
| `Halve`        | `MultiplyRate(0.5, ...)`       | —                 |
| `Double`       | `MultiplyRate(2.0, ...)`       | —                 |

Modifiers compose in declaration order.

## Range helpers

``` {.julia #subsidence-helpers}
_in_range(::Colon, _)                  = true
_in_range(r::Tuple{T,T}, v) where {T}  = r[1] <= v <= r[2]

function _t_overlap(m::AbstractSubsidenceModifier, t1::Time, t2::Time)
    m.t_range isa Colon && return (t1, t2)
    lo, hi = m.t_range
    a, b   = max(t1, lo), min(t2, hi)
    return a < b ? (a, b) : nothing
end
```

## Cumulative subsidence

The piecewise-constant modifier system is integrated analytically by walking
the union of event boundaries within $[t_0, t]$:

$$S(x, y, t) = \int_{t_0}^{t} r_{\mathrm{eff}}(x, y, \tau)\,\mathrm{d}\tau$$

``` {.julia #cumulative-subsidence}
"""
    cumulative_subsidence(base_rate_map, modifiers, x_axis, y_axis, t0)

Returns a closure `t -> Matrix{Length}` giving per-cell cumulative subsidence
from `t0` to `t`. Spatial masks are precomputed at construction time so
repeated evaluation is efficient.
"""
function cumulative_subsidence(
        base_rate_map::AbstractMatrix{<:Rate},
        modifiers::AbstractVector,
        x_axis::AbstractVector{<:Quantity},
        y_axis::AbstractVector{<:Quantity},
        t0::Time)

    nx, ny = size(base_rate_map)

    if isempty(modifiers)
        return (t::Time) -> base_rate_map .* (t - t0)
    end

    masks = [[_in_range(m.x_range, x_axis[i]) && _in_range(m.y_range, y_axis[j])
              for i in 1:nx, j in 1:ny]
             for m in modifiers]

    return function (t::Time)
        bounds = Time[t0, t]
        for m in modifiers
            ov = _t_overlap(m, t0, t)
            ov === nothing && continue
            push!(bounds, ov[1], ov[2])
        end
        unique!(sort!(bounds))

        accum = zeros(Length, nx, ny)
        for k in 1:length(bounds)-1
            ta, tb = bounds[k], bounds[k+1]
            tb <= ta && continue
            tmid = (ta + tb) / 2
            eff  = copy(base_rate_map)
            for (mi, m) in enumerate(modifiers)
                ov = _t_overlap(m, t0, t)
                ov === nothing && continue
                ov[1] <= tmid <= ov[2] || continue
                @inbounds for j in 1:ny, i in 1:nx
                    masks[mi][i, j] && (eff[i, j] = apply_rate(m, eff[i, j]))
                end
            end
            accum .+= eff .* (tb - ta)
        end
        return accum
    end
end
```

## Subsider

`subsider` applies one time step of subsidence to `state.bathymetry` — identical
pattern to the main branch (`bathymetry .-= Δσ`), but `Δσ` is now per-cell
for the matrix/modifier path:

``` {.julia #subsider}
function subsider(input::AbstractInput)
    # Scalar legacy path — bit-identical to main branch
    if input.subsidence_rate isa Quantity && isempty(input.subsidence_modifiers)
        Δσ = input.subsidence_rate * input.time.Δt
        return function (state::AbstractState)
            state.bathymetry .-= Δσ
        end
    end

    # Spatio-temporal path — same pattern, per-cell Δσ
    base     = subsidence_rate_map(input)
    x, y     = box_axes(input.box)
    Δt       = input.time.Δt
    t0       = input.time.t0
    get_time = time(input)
    cum      = cumulative_subsidence(base, input.subsidence_modifiers, x, y, t0)

    return function (state::AbstractState)
        t  = get_time(state)
        Δσ = cum(t + Δt) .- cum(t)
        state.bathymetry .-= Δσ
    end
end
```

## HDF5 serialization

The scalar `subsidence_rate` is always written for backward compatibility. The
full rate map and modifier descriptors are written only when non-uniform.

``` {.julia #waterdepth-write-header}
function write_header(input::AbstractInput, output::AbstractOutput)
    x, y = box_axes(input.box)
    t    = time_axis(input)
    set_attribute(output, "initial_topography", initial_topography(input) |> in_units_of(u"m"))
    set_attribute(output, "sea_level", input.sea_level.(t) .|> in_units_of(u"m"))

    rate_map   = subsidence_rate_map(input)
    scalar_rep = input.subsidence_rate isa Quantity ?
        input.subsidence_rate : sum(rate_map) / length(rate_map)
    set_attribute(output, "subsidence_rate", scalar_rep |> in_units_of(u"m/Myr"))

    if !(input.subsidence_rate isa Quantity)
        set_attribute(output, "subsidence_rate_map", rate_map .|> in_units_of(u"m/Myr"))
    end

    for (idx, m) in enumerate(input.subsidence_modifiers)
        prefix = "subsidence_modifiers/m$(idx)"
        set_attribute(output, "$(prefix)/kind", string(typeof(m).name.name))
        _write_modifier(output, prefix, m)
    end
end

_ser_axis(::Colon)  = [NaN, NaN]
_ser_axis(r::Tuple) = [ustrip(u"m",   r[1]), ustrip(u"m",   r[2])]
_ser_time(::Colon)  = [NaN, NaN]
_ser_time(r::Tuple) = [ustrip(u"Myr", r[1]), ustrip(u"Myr", r[2])]

function _write_modifier(out, prefix, m::MultiplyRate)
    set_attribute(out, "$(prefix)/factor",  m.factor)
    set_attribute(out, "$(prefix)/x_range", _ser_axis(m.x_range))
    set_attribute(out, "$(prefix)/y_range", _ser_axis(m.y_range))
    set_attribute(out, "$(prefix)/t_range", _ser_time(m.t_range))
end
function _write_modifier(out, prefix, m::AddRate)
    set_attribute(out, "$(prefix)/delta",   ustrip(u"m/Myr", m.delta))
    set_attribute(out, "$(prefix)/x_range", _ser_axis(m.x_range))
    set_attribute(out, "$(prefix)/y_range", _ser_axis(m.y_range))
    set_attribute(out, "$(prefix)/t_range", _ser_time(m.t_range))
end
function _write_modifier(out, prefix, m::SetRate)
    set_attribute(out, "$(prefix)/rate",    ustrip(u"m/Myr", m.rate))
    set_attribute(out, "$(prefix)/x_range", _ser_axis(m.x_range))
    set_attribute(out, "$(prefix)/y_range", _ser_axis(m.y_range))
    set_attribute(out, "$(prefix)/t_range", _ser_time(m.t_range))
end
```

## Deserialization

``` {.julia #subsidence-serialization}
deserialize_modifier(m::AbstractSubsidenceModifier) = m

function deserialize_modifier(d::AbstractDict)
    kind = d["kind"]
    xr   = _decode_axis(get(d, "x_range", [NaN, NaN]), u"m")
    yr   = _decode_axis(get(d, "y_range", [NaN, NaN]), u"m")
    tr   = _decode_axis(get(d, "t_range", [NaN, NaN]), u"Myr")
    if     kind == "MultiplyRate"
        return MultiplyRate(; factor=Float64(d["factor"]), x_range=xr, y_range=yr, t_range=tr)
    elseif kind == "AddRate"
        return AddRate(; delta=Float64(d["delta"])*u"m/Myr", x_range=xr, y_range=yr, t_range=tr)
    elseif kind == "SetRate"
        return SetRate(; rate=Float64(d["rate"])*u"m/Myr",   x_range=xr, y_range=yr, t_range=tr)
    else
        error("Subsidence: unknown modifier kind \"$(kind)\"")
    end
end

_decode_axis(v, unit) = (isnan(v[1]) || isnan(v[2])) ? (:) : (v[1]*unit, v[2]*unit)
```

## Example

Three ALCAP runs varying only the subsidence inputs:

``` {.julia .task file=examples/subsidence.jl}
#| creates:
#|   - data/output/subs-scalar.h5
#|   - data/output/subs-matrix.h5
#|   - data/output/subs-modifiers.h5
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
```

## Tests

``` {.julia file=test/Components/WaterDepthSpec.jl}
using CarboKitten
import CarboKitten.Components.WaterDepth as WD
using CarboKitten.Components.WaterDepth: MultiplyRate, Halve, cumulative_subsidence

# Original test from main — unchanged
@testset "Components/WaterDepth" begin
    input = WD.Input(
        box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
        time = TimeProperties(Δt=1.0u"Myr", steps=10),
        sea_level = t -> 2.0u"m",
        initial_topography = (x, y) -> -10.0u"m",
        subsidence_rate = 5.0u"m/Myr"
    )
    state = WD._initial_state(input)

    @test all(state.bathymetry .== WD.initial_topography(input))

    sub! = WD.subsider(input)
    sub!(state)
    @test all(state.bathymetry .== -15.0u"m")

    wd = WD.water_depth(input)
    @test all(wd(state) .== 17.0u"m")
end

@testset "Components/WaterDepth/matrix_rate" begin
    nx, ny = 5, 3
    rates  = [10.0u"m/Myr" + 10.0u"m/Myr" * (i - 1) for i in 1:nx, _ in 1:ny]
    input  = WD.Input(
        box = Box{Periodic{2}}(grid_size=(nx, ny), phys_scale=1000.0u"m"),
        time = TimeProperties(Δt=1.0u"Myr", steps=5),
        sea_level = t -> 0.0u"m",
        initial_topography = (x, y) -> -20.0u"m",
        subsidence_rate = rates)

    state = WD._initial_state(input)
    WD.subsider(input)(state)
    for i in 1:nx, j in 1:ny
        @test state.bathymetry[i, j] ≈ -20.0u"m" - rates[i, j] * 1.0u"Myr"
    end
end

@testset "Components/WaterDepth/modifier_halve" begin
    input = WD.Input(
        box = Box{Periodic{2}}(grid_size=(6, 2), phys_scale=1.0u"m"),
        time = TimeProperties(Δt=0.5u"Myr", steps=4),
        sea_level = t -> 0.0u"m",
        initial_topography = (x, y) -> -10.0u"m",
        subsidence_rate = 20.0u"m/Myr",
        subsidence_modifiers = [
            Halve(x_range=(0.0u"m", 3.0u"m"),
                  t_range=(0.0u"Myr", 0.5u"Myr"))])

    state = WD._initial_state(input)
    WD.subsider(input)(state)
    @test state.bathymetry[1, 1] ≈ -15.0u"m"   # halved: 10 m/Myr × 0.5 Myr
    @test state.bathymetry[6, 1] ≈ -20.0u"m"   # full:   20 m/Myr × 0.5 Myr
end

@testset "Components/WaterDepth/cumulative_analytic" begin
    base = fill(10.0u"m/Myr", 3, 3)
    x    = (1.0:3.0) * u"m"
    y    = (1.0:3.0) * u"m"
    cum  = cumulative_subsidence(base, [], x, y, 0.0u"Myr")
    @test all(cum(1.0u"Myr") .≈ 10.0u"m")
    @test all(cum(2.0u"Myr") .≈ 20.0u"m")

    # Halve for first 0.5 Myr: 5×0.5 + 10×0.5 = 7.5 m at t = 1 Myr
    mods = [MultiplyRate(0.5; t_range=(0.0u"Myr", 0.5u"Myr"))]
    cum2 = cumulative_subsidence(base, mods, x, y, 0.0u"Myr")
    @test all(cum2(1.0u"Myr") .≈ 7.5u"m")
end
```

## Component

``` {.julia file=src/Components/WaterDepth.jl}
@compose module WaterDepth
@mixin TimeIntegration, Boxes
using ..Common
using HDF5
using Unitful: ustrip
using ..TimeIntegration: time, time_axis
using CarboKitten.Components.Subsidence
import CarboKitten.Components.Subsidence: cumulative_subsidence

export water_depth, subsider, initial_topography, subsidence_rate_map
export AbstractSubsidenceModifier, MultiplyRate, AddRate, SetRate, Halve, Double
export cumulative_subsidence

<<waterdepth-input>>

<<subsidence-rate-map>>

<<subsider>>

@kwdef mutable struct State <: AbstractState
    bathymetry::Matrix{Height}
end

@constructor _initial_state(input)::State[bathymetry] =
    (bathymetry = initial_topography(input),)

function initial_state(input::AbstractInput)
    bathymetry = initial_topography(input)
    return State(step=0, bathymetry=bathymetry)
end

function initial_topography(input::AbstractInput)
    if input.initial_topography isa AbstractMatrix
        @assert size(input.initial_topography) == input.box.grid_size
        return input.initial_topography
    end

    x, y = box_axes(input.box)
    return input.initial_topography.(x, y')
end

function water_depth(input::AbstractInput)
    sea_level = input.sea_level
    get_time = time(input)

    return function (state::AbstractState)
        t = get_time(state)
        return sea_level(t) .- state.bathymetry
    end
end

<<waterdepth-write-header>>

end
```

``` {.julia file=src/Components/Subsidence.jl}
module Subsidence

using Unitful

export AbstractSubsidenceModifier, MultiplyRate, AddRate, SetRate, Halve, Double
export apply_rate, cumulative_subsidence, deserialize_modifier

const Rate     = typeof(1.0u"m/Myr")
const Time     = typeof(1.0u"Myr")
const Length   = typeof(1.0u"m")
const Location = typeof(1.0u"m")

const AxisSpec = Union{Colon,Tuple{Location,Location}}
const TimeSpec = Union{Colon,Tuple{Time,Time}}

abstract type AbstractSubsidenceModifier end

<<subsidence-modifiers>>

<<subsidence-helpers>>

<<cumulative-subsidence>>

<<subsidence-serialization>>

end
```
