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
