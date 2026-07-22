# ~/~ begin <<docs/src/components/waterdepth.md#src/Components/Subsidence.jl>>[init]
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

# ~/~ begin <<docs/src/components/waterdepth.md#subsidence-modifiers>>[init]
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
# ~/~ end

# ~/~ begin <<docs/src/components/waterdepth.md#subsidence-helpers>>[init]
_in_range(::Colon, _)                  = true
_in_range(r::Tuple{T,T}, v) where {T}  = r[1] <= v <= r[2]

function _t_overlap(m::AbstractSubsidenceModifier, t1::Time, t2::Time)
    m.t_range isa Colon && return (t1, t2)
    lo, hi = m.t_range
    a, b   = max(t1, lo), min(t2, hi)
    return a < b ? (a, b) : nothing
end
# ~/~ end

# ~/~ begin <<docs/src/components/waterdepth.md#cumulative-subsidence>>[init]
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
# ~/~ end

# ~/~ begin <<docs/src/components/waterdepth.md#subsidence-serialization>>[init]
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
# ~/~ end

end
# ~/~ end
