# ~/~ begin <<docs/src/facies-classification.md#examples/model/alcap/classify.jl>>[init]
#
# Post-deposition facies classification example
# ==============================================
# Reads the alcap-example.h5 produced by examples/model/alcap/run.jl,
# classifies the deposits into four lithological facies using Airy wave energy
# flux (W/m), water depth, and sediment proportion, then saves two fence
# diagrams:
#
#   data/output/alcap-example_fence_production.png   — production facies
#   data/output/alcap-example_fence_classified.png   — classified facies
#
# Run from the project root:
#   julia --project=. examples/model/alcap/classify.jl

module Script

using CairoMakie
using Makie
using Unitful
using CarboKitten
using CarboKitten.Export: read_volume
using CarboKitten.Visualization: fence_diagram!, fence_diagram
using CarboKitten.WaveField: AiryWaveField, AiryWaveComponent, energy_flux
using CarboKitten.FaciesClassification: FaciesRule, reclassify_data

const PATH = "data/output"
const TAG  = "alcap-example"

# ── Wave field ────────────────────────────────────────────────────────────────
# Use the same parameters as the ALCAP example run.
# Replace with your actual wave field if you used AiryWaveField in run.jl.
const WAVE_FIELD = AiryWaveField(components=[
    AiryWaveComponent(amplitude=1.5u"m", period=8.0u"s", direction=0.0),
])

# ── Energy thresholds (inspect before setting rules) ─────────────────────────
# energy_flux(WAVE_FIELD, d) for d in [5, 10, 20, 50] m:
#   ~3000 W/m, ~1500 W/m, ~600 W/m, ~200 W/m
const E_HIGH = energy_flux(WAVE_FIELD, 10.0u"m")   # ~1500 W/m: grainstone boundary
const E_MID  = energy_flux(WAVE_FIELD, 25.0u"m")   # ~400 W/m: packstone boundary

# ── Classification rules ──────────────────────────────────────────────────────
# Production facies: 1 = euphotic, 2 = oligophotic, 3 = aphotic

const RULES = [
    # Grainstone: shallow, wave-agitated, euphotic-dominated
    FaciesRule(
        name               = "grainstone",
        sediment_fractions = Dict(1 => (0.5, 1.0)),
        depth_range        = (0.0u"m",  20.0u"m"),
        wave_energy_range  = (E_HIGH, Inf*u"W/m")),

    # Packstone: moderate energy, euphotic/oligophotic mix
    FaciesRule(
        name               = "packstone",
        sediment_fractions = Dict(1 => (0.2, 0.8)),
        depth_range        = (0.0u"m",  40.0u"m"),
        wave_energy_range  = (E_MID, Inf*u"W/m")),

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

const RULE_LABELS = [[r.name for r in RULES]; "unclassified"]

# ── Fence slices ─────────────────────────────────────────────────────────────
const X_SLICES = [10, 30, 50]
const Y_SLICES = [2.0u"km", 4.0u"km", 6.0u"km"]

function _legend!(fig, pos, n, labels, title)
    colors   = Makie.wong_colors()[1:n]
    elements = [PolyElement(color=colors[i]) for i in 1:n]
    Legend(fig[pos...], elements, labels, title)
end

function _fence_fig(header, slices, title, n_facies, labels)
    fig = Figure(size=(1400, 800))
    ax  = Axis3(fig[1, 1]; title=title)
    ax.azimuth   = -π/3
    ax.elevation =  π/8
    fence_diagram!(ax, header, slices;
        show_unconformities = 10,
        show_coeval_lines   = true,
        show_sealevel       = true,
        color_by            = :facies)
    _legend!(fig, (1, 2), n_facies, labels, "Facies")
    return fig
end

function main()
    filename = "$(PATH)/$(TAG).h5"
    @info "Reading $filename …"
    header, vol = read_volume(filename, :topography)
    nx, ny = header.grid_size

    # Build slice list (same geometry for both panels)
    x_idx = [i for i in X_SLICES if i <= nx]
    y_idx = [argmin(abs.(header.axes.y .- p))
              for p in Y_SLICES if uconvert(unit(header.axes.y[1]), p) <= header.axes.y[end]]
    slices_prod = vcat([vol[i, :] for i in x_idx],
                       [vol[:, j] for j in y_idx])
    isempty(slices_prod) && error("No slices selected — check X_SLICES/Y_SLICES against grid size $(nx)×$(ny)")
    @info "$(length(slices_prod)) slices selected, types: $(typeof.(slices_prod))"

    # Classify each slice
    @info "Classifying …"
    classified = [reclassify_data(header, s, RULES; wave_field=WAVE_FIELD)
                  for s in slices_prod]
    new_header  = classified[1][1]
    slices_cls  = [c[2] for c in classified]

    prod_labels = ["euphotic", "oligophotic", "aphotic"]   # from run.jl

    # Production-facies fence diagram
    fig_prod = _fence_fig(header, slices_prod,
                          "Production facies", header.n_facies, prod_labels)
    out_prod = "$(PATH)/$(TAG)_fence_production.png"
    save(out_prod, fig_prod)
    @info "Saved $out_prod"

    # Classified-facies fence diagram
    fig_cls = _fence_fig(new_header, slices_cls,
                         "Classified facies", new_header.n_facies, RULE_LABELS)
    out_cls = "$(PATH)/$(TAG)_fence_classified.png"
    save(out_cls, fig_cls)
    @info "Saved $out_cls"
end

end  # module Script

Script.main()
# ~/~ end
