# ~/~ begin <<docs/src/facies-classification.md#src/FaciesClassification.jl>>[init]
# ~/~ begin <<docs/src/facies-classification.md#src/FaciesClassification.jl>>[init]
"""
    FaciesClassification
...
"""
module FaciesClassification

using Unitful
using ..Output.Abstract: Data, Header, water_depth
using ..WaveField: AiryWaveField, energy_flux

export FaciesRule, classify_block, reclassify_data

# ~/~ begin <<docs/src/facies-classification.md#facies-rule>>[init]
"""
    FaciesRule(; name, sediment_fractions, depth_range, wave_energy_range)

A single rule for post-depositional facies classification.

## Fields
- `name::String` — label assigned to blocks that match this rule.
- `sediment_fractions::Union{Nothing, Dict{Int,NTuple{2,Float64}}}` —
  optional per-production-facies fraction constraints.  Keys are 1-based
  production-facies indices; values are `(lo, hi)` with `0 ≤ lo ≤ hi ≤ 1`.
  All listed constraints must be satisfied simultaneously.  `nothing` means
  no constraint on sediment proportion.
- `depth_range::NTuple{2,<:Quantity}` — `(min_depth, max_depth)` water-depth
  window in metres.  Bounds are inclusive.  Use `(-Inf*u"m", Inf*u"m")` to
  leave depth unconstrained.
- `wave_energy_range::NTuple{2,<:Quantity}` — `(lo, hi)` for Airy wave energy
  flux in W/m, from `WaveField.energy_flux`.  Use `(0.0u"W/m", Inf*u"W/m")`
  to leave wave energy unconstrained.

Rules are evaluated in order; first match wins.  Unmatched blocks go to the
fallback class at index `length(rules) + 1`.
"""
@kwdef struct FaciesRule
    name::String
    sediment_fractions::Union{Nothing,Dict{Int,NTuple{2,Float64}}} = nothing
    depth_range::NTuple{2,<:Quantity}       = (-Inf * u"m",   Inf * u"m")
    wave_energy_range::NTuple{2,<:Quantity} = (0.0  * u"W/m", Inf * u"W/m")
end
# ~/~ end

# ~/~ begin <<docs/src/facies-classification.md#classify-block>>[init]
"""
    classify_block(rules, fractions, depth, wave_energy) -> Int

Classify a single deposited block against an ordered list of `FaciesRule`s.

Returns the 1-based index of the first matching rule, or `length(rules) + 1`
(fallback) when no rule matches.
"""
function classify_block(rules::AbstractVector{FaciesRule},
                        fractions::AbstractVector{<:Real},
                        depth::Quantity,
                        wave_energy::Quantity)::Int
    for (i, rule) in enumerate(rules)
        depth < rule.depth_range[1] && continue
        depth > rule.depth_range[2] && continue
        wave_energy < rule.wave_energy_range[1] && continue
        wave_energy > rule.wave_energy_range[2] && continue
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
# ~/~ end

# ~/~ begin <<docs/src/facies-classification.md#reclassify-data>>[init]
"""
    reclassify_data(header, data, rules; wave_field=nothing) -> (Header, Data)

Reclassify a CarboKitten `Data` object using an ordered vector of `FaciesRule`s.

`wave_field::Union{AiryWaveField, Nothing}` — the same wave field used during
the model run.  `WaveField.energy_flux(wave_field, wd)` is called at each cell
to obtain the wave energy in W/m.  When `nothing`, wave energy is `0 W/m`
everywhere.

Water depth is read from `data.water_depth` when stored (run with
`save_water_depth=true`), otherwise reconstructed from the header.

Returns `(new_header, new_data)` with the same `Data{F,D}` type, fully
compatible with all existing visualisation routines.
"""
function reclassify_data(header::Header,
                         data::Data{F,D},
                         rules::AbstractVector{FaciesRule};
                         wave_field::Union{AiryWaveField,Nothing}=nothing) where {F,D}
    n_prod  = size(data.deposition, 1)    # original production-facies count
    n_class = length(rules) + 1           # classified facies + fallback
    dep_sz  = size(data.deposition)
    sp_size = dep_sz[2:end-1]             # (): column, (nx,): slice, (nx,ny): volume
    n_t     = dep_sz[end]

    # --- allocate output arrays (same element type / units as input) ---
    zero_like(a) = zeros(eltype(a), n_class, sp_size..., n_t)
    dep_out  = zero_like(data.deposition)
    prod_out = zero_like(data.production)
    dis_out  = zero_like(data.disintegration)

    # --- water depth: stored field preferred, fallback to header reconstruction ---
    wd_array = water_depth(header, data)   # shape (spatial..., n_t)

    _wd(sp_idx::CartesianIndex{0}, t) = wd_array[t]
    _wd(sp_idx, t)                    = wd_array[sp_idx, t]

    zero_energy = 0.0u"W/m"

    # --- iterate over every spatial cell and time step ---
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

            we  = wave_field !== nothing ? energy_flux(wave_field, wd_val) : zero_energy
            cls = classify_block(rules, fractions, wd_val, we)

            dep_out[cls,  sp_idx, t_idx] += total_dep
            prod_out[cls, sp_idx, t_idx] += sum(prod_col)
            dis_out[cls,  sp_idx, t_idx] += sum(dis_col)
        end
    end

    # --- build new header ---
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
        active_layer       = nothing,
        water_depth        = data.water_depth)

    return new_header, new_data
end
# ~/~ end

end  # module FaciesClassification
# ~/~ end
# ~/~ end
