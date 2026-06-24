# Post-Deposition Facies Classification

CarboKitten models record deposition in terms of **production facies** —
categories that reflect the biological or abiotic source of the sediment
(e.g. euphotic, oligophotic, aphotic carbonate).  In practice, stratigraphers
classify sediment packages into **depositional facies** (wackestone, packstone,
grainstone, …) based on observable proxies: sediment composition, palaeo-water
depth, and wave energy.

The `FaciesClassification` module provides a lightweight, non-invasive workflow
that performs this second classification step entirely in post-processing, on
the `Data` objects returned by the standard read routines.  No model, component,
or output writer is modified.  The approach is conceptually analogous to the
rule-based facies assignment in DIONISOSFlow.

Wave energy is computed using the physically consistent Airy wave theory
formulation provided by the `WaveField` module (see [Wave Field](wavefield.md)).
The same `AiryWaveField` used to drive sediment transport during the run is
passed to `reclassify_data` at classification time — there is no separate proxy
to calibrate.

## Design overview

```
                        ┌──────────────────────────────┐
  run_model  ──────────▶│  Data  (production facies)    │
                        │  data.water_depth  (optional) │
                        └────────────┬─────────────────┘
                                     │ reclassify_data(rules; wave_field)
                        ┌────────────▼─────────────────┐
                        │  Data  (classified facies)    │──▶ fence_diagram
                        └──────────────────────────────┘     sediment_profile …
```

The input `Data` object carries `n_prod` production facies; the output carries
`n_class = length(rules) + 1` classified facies (the last slot is always the
fallback for blocks matching no rule).  Because the return type is the same
`Data{F,D}` struct, **all existing visualisation routines work without
modification**.

### Water depth

`reclassify_data` obtains water depth through
`Output.Abstract.water_depth(header, data)`, which has two paths:

1. **Stored field** (`data.water_depth !== nothing`): returned directly.
   Required when subsidence is spatially variable.  Enable with
   `save_water_depth = true` in the model input.
2. **Reconstructed from header** (fallback): `sl(t) - h₀ - Δh(t) + D·t`,
   exact for spatially uniform subsidence.

### Wave energy

Wave energy flux at each cell is computed as `WaveField.energy_flux(wave_field,
wd)` in W/m.  This is the same quantity as in the `WaveField` module —
`(1/8) ρ g H_eff² c_g` summed over all wave components — with depth-limited
breaking applied.  If no wave field is supplied, wave energy is `0 W/m`
everywhere and `wave_energy_range` rules never fire.

## FaciesRule

``` {.julia #facies-rule}
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
```

## classify_block

``` {.julia #classify-block}
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
```

## reclassify_data

``` {.julia #reclassify-data}
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
    # ... (full implementation in src/FaciesClassification.jl)
end
```

## User-facing example

```julia
using CarboKitten
using CarboKitten.Export: read_volume
using CarboKitten.WaveField: AiryWaveField, AiryWaveComponent, energy_flux
using CarboKitten.FaciesClassification: FaciesRule, reclassify_data
using Unitful

# ── 1. Define wave field (same as used in the model run) ─────────────────────
wf = AiryWaveField(components=[
    AiryWaveComponent(amplitude=1.5u"m", period=8.0u"s", direction=0.0),
    AiryWaveComponent(amplitude=0.5u"m", period=5.0u"s", direction=π/4),
])

# Inspect thresholds: energy at candidate depth boundaries
for d in [5.0, 10.0, 20.0, 50.0]
    @info "E at $(d)m = $(energy_flux(wf, d*u\"m\"))"
end

# ── 2. Run model (with water depth stored for variable-subsidence safety) ────
input = ALCAP.Input(
    # ... existing fields ...
    save_water_depth = true,
    facies = [ALCAP.Facies(..., wave_velocity=wf), ...])

run_model(Model{ALCAP}, input, "output/run.h5")

# ── 3. Read output ────────────────────────────────────────────────────────────
header, vol = read_volume("output/run.h5", :topography)

# ── 4. Define classification rules ───────────────────────────────────────────
# Production facies: 1 = euphotic, 2 = oligophotic, 3 = aphotic
# Energy thresholds chosen from the inspection above.

rules = [
    # Grainstone: shallow, high wave energy, euphotic-dominated
    FaciesRule(
        name               = "grainstone",
        sediment_fractions = Dict(1 => (0.5, 1.0)),
        depth_range        = (0.0u"m",  20.0u"m"),
        wave_energy_range  = (500.0u"W/m", Inf*u"W/m")),

    # Packstone: moderate energy, mixed euphotic/oligophotic
    FaciesRule(
        name               = "packstone",
        sediment_fractions = Dict(1 => (0.2, 0.8)),
        depth_range        = (0.0u"m",  40.0u"m")),

    # Wackestone: low energy, oligophotic-dominated
    FaciesRule(
        name               = "wackestone",
        sediment_fractions = Dict(2 => (0.3, 1.0)),
        depth_range        = (5.0u"m",  80.0u"m")),

    # Mudstone: deep, aphotic-dominated
    FaciesRule(
        name        = "mudstone",
        depth_range = (10.0u"m", Inf*u"m")),
    # Index 5 = fallback (unclassified)
]

# ── 5. Classify ───────────────────────────────────────────────────────────────
new_header, new_vol = reclassify_data(header, vol, rules; wave_field=wf)

# new_header.attributes["classified_facies"]
# => ["grainstone", "packstone", "wackestone", "mudstone", "fallback"]

# ── 6. Visualise (unchanged routines) ────────────────────────────────────────
using CairoMakie, CarboKitten.Visualization
fig = fence_diagram(new_header, new_vol;
    x_slices=[10, 30, 50], y_slices=[2.0u"km", 4.0u"km"])
save("fence_classified.png", fig)
```

## Setting wave energy thresholds

Unlike depth, wave energy thresholds are not immediately intuitive.  The
suggested workflow is to call `energy_flux(wf, d*u"m")` for a range of depths
before setting rules, as shown in the example above.  For the default ALCAP
example wave field (1.5 m / 8 s swell):

| Depth | Energy flux |
|-------|-------------|
| 5 m   | ~3 000 W/m  |
| 10 m  | ~1 500 W/m  |
| 20 m  | ~600 W/m    |
| 50 m  | ~200 W/m    |

These values scale with `H²`, so a 3 m wave has roughly 4× the flux of a 1.5 m
wave at the same depth.

## Limitations

- All production at one cell in one time step is collapsed to a single
  classified-facies bucket.  Intra-timestep mixing is not resolved.
- `active_layer` is not reclassified (set to `nothing`): it represents
  pre-lithification mobile sediment, not a deposited block.
- Wave energy is evaluated at the *palaeo-water depth* of the deposited block,
  not at present-day depth.  This is the intended behaviour.
- When `save_water_depth = false` and subsidence is spatially variable,
  classification results will be incorrect.

## Source file

``` {.julia file=src/FaciesClassification.jl}
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

<<facies-rule>>

<<classify-block>>

<<reclassify-data>>

end  # module FaciesClassification
# ~/~ end
```

## Test specification

``` {.julia file=test/FaciesClassificationSpec.jl}
# ~/~ begin <<docs/src/facies-classification.md#test/FaciesClassificationSpec.jl>>[init]
# ... (see test/FaciesClassificationSpec.jl)
# ~/~ end
```
