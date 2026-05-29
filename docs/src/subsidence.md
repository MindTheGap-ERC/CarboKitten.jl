# Subsidence

Subsidence modifiers and the cumulative-subsidence integral.
 
This module is intentionally low-level — it depends only on `Unitful` — so
both the runtime path (`Components.WaterDepth`) and the post-hoc analysis
path (`Output.Abstract`, plot routines) can share the same algorithm.
 
## Modifiers
 
A modifier locally alters the subsidence rate inside an (x, y, t) box:
 
```julia
MultiplyRate(0.5;  x_range=(0u"m", 1u"km"), t_range=(0u"Myr", 0.5u"Myr"))
AddRate(20u"m/Myr"; x_range=(3u"km", 6u"km"))
SetRate(0u"m/Myr";  x_range=(12u"km", 15u"km"), t_range=(0.75u"Myr", 1u"Myr"))
Halve(t_range=(0u"Myr", 0.2u"Myr"))     # MultiplyRate(0.5)
Double(t_range=(0u"Myr", 0.2u"Myr"))    # MultiplyRate(2.0)
```
 
`x_range`, `y_range`, `t_range` accept `:` for "no restriction", or a 2-tuple
of `Quantity` interpreted as an inclusive `[lo, hi]` interval.
 
## Cumulative-subsidence integral
 
For a base rate map (Matrix{Rate}) and a list of modifiers, the cumulative
subsidence at time `t` is the piecewise integral
 
```math
S(x, y, t) = ∫_{t_0}^{t} r_eff(x, y, τ) dτ
```
 
where `r_eff` is `base_rate` sequentially transformed by every modifier whose
(x, y, t)-box contains the point. Modifiers apply in declaration order, each
operating on the rate produced by the previous one.

## Compatibility

 To ensure compatibility with the plotting routines, these have been updated to receive cumulative subsidence. This allows the plots to account for time- and space-dependent subsidence when calculating the water depth and stratigraphic height.

### Tests - Subsidence-modification examples

The following examples consist of three runs of the ALCAP example, each varying only the subsidence inputs.
All three should work side-by-side and produce H5 outputs.

1. Legacy: runs with a uniform scalar subsidence rate (50 m/Myr everywhere).
2. Per-cell rate: runs with a ramped subsidence increasing along x.
3. With modifiers: runs with a uniform base rate, which is halved for a specific region over a specific time window while another region receives an additive bump.

```{.julia .task file=examples/subsidence.jl}
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
                x_range=(3000.0u"m", 6000.0u"m"),
                y_range=(2000.0u"m", 5000.0u"m")),
        # Pin the rate to 0 over the rightmost strip for the last quarter
        SetRate(0.0u"m/Myr";
                x_range=(12000.0u"m", 15000.0u"m"),
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
```

### Implementation

```{.julia file=src/Subsidence.jl}
module Subsidence

using Unitful

export AbstractSubsidenceModifier, MultiplyRate, AddRate, SetRate, Halve, Double
export apply_rate, cumulative_subsidence, deserialize_modifier

const Rate     = typeof(1.0u"m/Myr")
const Time     = typeof(1.0u"Myr")
const Length   = typeof(1.0u"m")
const Location = typeof(1.0u"m")

# Specs for the (x, y) and t axes of a modifier's box.
const AxisSpec = Union{Colon,Tuple{Location,Location},AbstractRange{<:Location}}
const TimeSpec = Union{Colon,Tuple{Time,Time},AbstractRange{<:Time}}

abstract type AbstractSubsidenceModifier end

"""
    MultiplyRate(factor; x_range=:, y_range=:, t_range=:)

Multiply the subsidence rate by `factor` inside the given box.
"""
@kwdef struct MultiplyRate <: AbstractSubsidenceModifier
    factor::Float64
    x_range::AxisSpec = (:)
    y_range::AxisSpec = (:)
    t_range::TimeSpec = (:)
end
MultiplyRate(factor::Real; kwargs...) = MultiplyRate(; factor=Float64(factor), kwargs...)

Halve(;  kwargs...) = MultiplyRate(0.5; kwargs...)
Double(; kwargs...) = MultiplyRate(2.0; kwargs...)

"""
    AddRate(delta; x_range=:, y_range=:, t_range=:)

Add `delta` (a `Rate`) to the subsidence rate inside the given box. Use a
negative `delta` for subtraction.
"""
@kwdef struct AddRate <: AbstractSubsidenceModifier
    delta::Rate
    x_range::AxisSpec = (:)
    y_range::AxisSpec = (:)
    t_range::TimeSpec = (:)
end
AddRate(delta::Rate; kwargs...) = AddRate(; delta=delta, kwargs...)

"""
    SetRate(rate; x_range=:, y_range=:, t_range=:)

Override the subsidence rate to `rate` (a `Rate`) inside the given box.
"""
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

# -- range membership ---------------------------------------------------------
_in_range(::Colon, _)                 = true
_in_range(r::AbstractRange, v)        = first(r) <= v <= last(r)
_in_range(r::Tuple{T,T}, v) where {T} = r[1] <= v <= r[2]

"""
    modifier_active(m, x, y, t) -> Bool

Whether modifier `m` applies at physical location `(x, y)` and time `t`.
"""
modifier_active(m::AbstractSubsidenceModifier, x, y, t) =
    _in_range(m.x_range, x) && _in_range(m.y_range, y) && _in_range(m.t_range, t)

# Project the modifier's t_range onto [t1, t2]; nothing if no overlap.
function _t_overlap(m::AbstractSubsidenceModifier, t1::Time, t2::Time)
    if m.t_range isa Colon
        return (t1, t2)
    elseif m.t_range isa Tuple
        lo, hi = m.t_range
    else
        lo, hi = first(m.t_range), last(m.t_range)
    end
    a = max(t1, lo)
    b = min(t2, hi)
    return a < b ? (a, b) : nothing
end

# =============================================================================
# Cumulative-subsidence integral
# =============================================================================

"""
    cumulative_subsidence(base_rate_map, modifiers, x_axis, y_axis, t0) -> closure

Build a closure `t -> Matrix{Length}` that returns the cumulative subsidence
from `t0` to `t`, per cell, over the grid defined by `x_axis × y_axis`. The
closure precomputes per-modifier spatial masks so it is fast to call many
times (e.g. once per write step in a plot routine).

    cumulative_subsidence(base_rate_map, modifiers, x_axis, y_axis, t0, t) -> Matrix{Length}

Convenience form that evaluates immediately.
"""
function cumulative_subsidence(
        base_rate_map::AbstractMatrix{<:Rate},
        modifiers::AbstractVector,
        x_axis::AbstractVector{<:Quantity},
        y_axis::AbstractVector{<:Quantity},
        t0::Time)
    nx, ny = size(base_rate_map)
    @assert length(x_axis) == nx "x_axis length $(length(x_axis)) does not match base_rate_map first dim $(nx)"
    @assert length(y_axis) == ny "y_axis length $(length(y_axis)) does not match base_rate_map second dim $(ny)"

    # Fast path: no modifiers
    if isempty(modifiers)
        return function (t::Time)
            return base_rate_map .* (t - t0)
        end
    end

    # Precompute static spatial masks for each modifier
    spatial_masks = [
        [_in_range(m.x_range, x_axis[i]) && _in_range(m.y_range, y_axis[j])
         for i in 1:nx, j in 1:ny]
        for m in modifiers]

    return function (t::Time)
        # Walk modifier event boundaries within [t0, t] for piecewise integration.
        boundaries = Time[t0, t]
        for m in modifiers
            ov = _t_overlap(m, t0, t)
            ov === nothing && continue
            push!(boundaries, ov[1]); push!(boundaries, ov[2])
        end
        sort!(boundaries)
        unique!(boundaries)

        accum = zeros(Length, nx, ny)
        for k in 1:length(boundaries)-1
            ta, tb = boundaries[k], boundaries[k+1]
            tb <= ta && continue
            tmid = (ta + tb) / 2

            # Effective rate on the grid for this sub-interval: base rate, then
            # each active modifier applied in order.
            eff = copy(base_rate_map)
            for (mi, m) in enumerate(modifiers)
                ov = _t_overlap(m, t0, t)
                ov === nothing && continue
                if ov[1] <= tmid <= ov[2]
                    mask = spatial_masks[mi]
                    @inbounds for j in 1:ny, i in 1:nx
                        if mask[i, j]
                            eff[i, j] = apply_rate(m, eff[i, j])
                        end
                    end
                end
            end
            accum .+= eff .* (tb - ta)
        end
        return accum
    end
end

cumulative_subsidence(base, modifiers, x_axis, y_axis, t0::Time, t::Time) =
    cumulative_subsidence(base, modifiers, x_axis, y_axis, t0)(t)

# =============================================================================
# (De)serialization
# =============================================================================

# Already-live modifier — identity.
deserialize_modifier(m::AbstractSubsidenceModifier) = m

"""
    deserialize_modifier(d)

Reconstruct a modifier object from a dict-like descriptor (as produced when
reading from HDF5). The descriptor must have a `"kind"` key naming the
concrete type ("MultiplyRate", "AddRate", "SetRate") and the appropriate
parameter keys.

`x_range`, `y_range`, `t_range` are stored as 2-element `Float64` arrays; a
pair of `NaN`s decodes to `:` (no restriction).
"""
function deserialize_modifier(d::AbstractDict)
    kind = d["kind"]
    xr = _decode_axis(get(d, "x_range", [NaN, NaN]), u"m")
    yr = _decode_axis(get(d, "y_range", [NaN, NaN]), u"m")
    tr = _decode_axis(get(d, "t_range", [NaN, NaN]), u"Myr")
    if kind == "MultiplyRate"
        return MultiplyRate(factor=Float64(d["factor"]), x_range=xr, y_range=yr, t_range=tr)
    elseif kind == "AddRate"
        return AddRate(delta=Float64(d["delta"])*u"m/Myr", x_range=xr, y_range=yr, t_range=tr)
    elseif kind == "SetRate"
        return SetRate(rate=Float64(d["rate"])*u"m/Myr", x_range=xr, y_range=yr, t_range=tr)
    else
        error("Subsidence: unknown modifier kind \"$(kind)\"")
    end
end

function _decode_axis(v::AbstractVector, unit)
    a, b = v[1], v[2]
    if isnan(a) || isnan(b)
        return (:)
    end
    return (a * unit, b * unit)
end

end



```
