# ~/~ begin <<docs/src/subsidence.md#src/Subsidence.jl>>[init]
module Subsidence

using Unitful

export AbstractSubsidenceModifier, MultiplyRate, AddRate, SetRate, Halve, Double
export apply_rate, cumulative_subsidence, deserialize_modifier

const Rate     = typeof(1.0u"m/Myr")
const Time     = typeof(1.0u"Myr")
const Length   = typeof(1.0u"m")
const Location = typeof(1.0u"m")

#Normalize  constants to expected units.
_rate(r::Quantity) = uconvert(u"m/Myr", r)
_length(x::Quantity) = uconvert(u"m", x)
_time(t::Quantity) = uconvert(u"Myr", t)

_normalize_axis(::Colon) = (:)
_normalize_axis(r::Tuple{<:Quantity,<:Quantity}) =
    (_length(r[1]), _length(r[2]))
_normalize_axis(r::AbstractRange{<:Quantity}) =
    _length.(r)

_normalize_time(::Colon) = (:)
_normalize_time(r::Tuple{<:Quantity,<:Quantity}) =
    (_time(r[1]), _time(r[2]))
_normalize_time(r::AbstractRange{<:Quantity}) =
    _time.(r)

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
MultiplyRate(factor::Real; x_range = (:), y_range = (:), t_range = (:)) = MultiplyRate(factor = Float64(factor), x_range = _normalize_axis(x_range), y_range = _normalize_axis(y_range), t_range = _normalize_time(t_range))

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
AddRate(delta::Quantity; x_range = (:), y_range = (:), t_range = (:)) = AddRate(delta = _rate(delta), x_range = _normalize_axis(x_range), y_range = _normalize_axis(y_range), t_range = _normalize_time(t_range))

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
SetRate(rate::Quantity; x_range = (:), y_range = (:), t_range = (:)) = SetRate( rate = _rate(rate),  x_range = _normalize_axis(x_range), y_range = _normalize_axis(y_range), t_range = _normalize_time(t_range))

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



# ~/~ end
