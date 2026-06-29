# run-ecl.jl — Upper Malm, Eclépens
#
# Usage:
#   julia --project=. examples/model/run-ecl.jl

module Script

using CarboKitten
using CarboKitten.Models: WithoutCA
using CarboKitten.Production: InterpolatedProduction, MultiplyProduction
using CarboKitten.WaveField: AiryWaveField, AiryWaveComponent
using CarboKitten.FaciesClassification: FaciesRule, reclassify_data,
                                         reclassify_volume, save_classified
using CarboKitten.Export: read_volume, read_slice
using CarboKitten.Visualization: fence_diagram, fence_diagram!
using CSV, DataFrames
using Interpolations
using CairoMakie
using Unitful
using Random, Statistics

include(joinpath(@__DIR__, "../../ext/MapView.jl"))
include(joinpath(@__DIR__, "../../ext/StratigraphicColumn.jl"))
using .MapView: map_view
using .StratigraphicColumn: stratigraphic_column!

# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────

const BATHY_CSV       = "data/input/bathymetry.csv"
const SEALEVEL_CSV    = "data/input/sealevel.csv"
const OUTPUT_FILE     = "data/output/eclepens_withoutca.h5"
const OUTPUT_CLS_FILE = "data/output/eclepens_withoutca_classified.h5"
const FIG_DIR         = "data/output/eclepens_withoutca_figs"

# ─────────────────────────────────────────────────────────────────────────────
# 1. Bathymetry — reduce ooid blankets, add upper-platform coral micro-relief
# ─────────────────────────────────────────────────────────────────────────────

bathy_raw = Matrix{Float64}(CSV.read(BATHY_CSV, DataFrame; header=false))
const RAW_GRID_SIZE = (size(bathy_raw, 1), size(bathy_raw, 2))

const BATHY_TARGET_MIN =  8.8
const BATHY_TARGET_MAX = 66.0
const BATHY_GAMMA      = 1.35

function gaussian_relief(grid_size;
                         n::Int,
                         amp_range::Tuple{Float64, Float64},
                         rx_range::Tuple{Float64, Float64},
                         ry_range::Tuple{Float64, Float64},
                         seed::Int)
    rng = MersenneTwister(seed)
    nx, ny = grid_size
    relief = zeros(Float64, nx, ny)

    for _ in 1:n
        cx  = 1.0 + (nx - 1.0) * rand(rng)
        cy  = 1.0 + (ny - 1.0) * rand(rng)
        amp = amp_range[1] + (amp_range[2] - amp_range[1]) * rand(rng)
        rx  = rx_range[1]  + (rx_range[2]  - rx_range[1])  * rand(rng)
        ry  = ry_range[1]  + (ry_range[2]  - ry_range[1])  * rand(rng)

        for j in 1:ny, i in 1:nx
            dx  = i - cx
            dy0 = abs(j - cy)
            dy  = min(dy0, ny - dy0)
            relief[i, j] += amp * exp(-0.5 * ((dx / rx)^2 + (dy / ry)^2))
        end
    end

    return relief
end

let qlo = quantile(vec(bathy_raw), 0.01),
    qhi = quantile(vec(bathy_raw), 0.99)

    z = clamp.((bathy_raw .- qlo) ./ (qhi - qlo), 0.0, 1.0)

    bathy_base = BATHY_TARGET_MIN .+
                 (BATHY_TARGET_MAX - BATHY_TARGET_MIN) .* z .^ BATHY_GAMMA

    # Fewer and weaker broad highs: ooids stay as shoal patches, not sheets.
        large_shoal_highs = gaussian_relief(
        RAW_GRID_SIZE;
        n         = max(22, round(Int, prod(RAW_GRID_SIZE) / 3000)),
        amp_range = (3.0, 6.5),
        rx_range  = (1.6, 3.6),
        ry_range  = (1.8, 5.0),
        seed      = 1573,
    )

    # More small upper-platform knobs: creates coral opportunities on the
    # upper half of the platform without using CA or masks.
        coral_knobs = gaussian_relief(
        RAW_GRID_SIZE;
        n         = max(80, round(Int, prod(RAW_GRID_SIZE) / 850)),
        amp_range = (0.4, 1.6),
        rx_range  = (0.7, 1.5),
        ry_range  = (0.7, 1.8),
        seed      = 1450,
    )

    # Broad shallow pans preserve restricted lagoon around shoals.
    lagoon_pans = gaussian_relief(
        RAW_GRID_SIZE;
        n         = max(12, round(Int, prod(RAW_GRID_SIZE) / 6000)),
        amp_range = (1.0, 2.8),
        rx_range  = (8.0, 18.0),
        ry_range  = (8.0, 18.0),
        seed      = 1512,
    )

    # Small bathymetric lows split otherwise connected shallow ooid areas.
    shoal_breaks = gaussian_relief(
        RAW_GRID_SIZE;
        n         = max(45, round(Int, prod(RAW_GRID_SIZE) / 1300)),
        amp_range = (0.7, 2.0),
        rx_range  = (3.5, 9.0),
        ry_range  = (0.6, 1.8),
        seed      = 1601,
    )

    coral_knob_gate = clamp.((bathy_base .- 10.0) ./ 12.0, 0.0, 1.0)

    global bathy_rescaled = clamp.(
        bathy_base
        .- large_shoal_highs
        .- coral_knob_gate .* coral_knobs
        .+ lagoon_pans
        .+ shoal_breaks,
        2.8,
        BATHY_TARGET_MAX,
    )
    end

bathy_m      = bathy_rescaled * u"m"
initial_topo = -bathy_m
const GRID_SIZE = (size(bathy_m, 1), size(bathy_m, 2))

@info "Bathymetry: $(GRID_SIZE[1])×$(GRID_SIZE[2]), textured depth $(minimum(bathy_rescaled))–$(maximum(bathy_rescaled)) m"
# ─────────────────────────────────────────────────────────────────────────────
# 2. Sea-level curve — unchanged
# ─────────────────────────────────────────────────────────────────────────────

const T_START_MA = -157.3
const T_END_MA   = -145.0

sl_df  = CSV.read(SEALEVEL_CSV, DataFrame)
age_Ma = Float64.(sl_df.age_Ma)
sl_z   = Float64.(sl_df.sealevel_m)

@assert T_START_MA >= minimum(age_Ma)
@assert T_END_MA   <= maximum(age_Ma)
@assert T_END_MA   >  T_START_MA

t_elapsed = age_Ma .- T_START_MA
sl_itp    = linear_interpolation(t_elapsed, sl_z; extrapolation_bc=Flat())
sl_at_t0  = sl_itp(0.0)
sea_level = t -> (sl_itp(ustrip(u"Myr", t)) - sl_at_t0) * u"m"

# ─────────────────────────────────────────────────────────────────────────────
# 3. Wave field — classification only
# ─────────────────────────────────────────────────────────────────────────────

wave_field = AiryWaveField(components=[
    AiryWaveComponent(amplitude=2.225u"m", period=12.12u"s", direction=deg2rad(245)),
    AiryWaveComponent(amplitude=1.0u"m",   period=10.0u"s",  direction=deg2rad(60)),
])

# ─────────────────────────────────────────────────────────────────────────────
# 4. Production facies
#
# RATES: Dionisos literature values (Sultana et al. 2022):
#   ooid=50, coral=10, mud=8, grains=12, sabkha≈0 m/Myr
#
# PRODUCTION BUDGET vs accommodation (≈20 m/Myr):
#
#   Zone         Depth    Budget (m/Myr)   Behaviour
#   ──────────────────────────────────────────────────────────────
#   Ooid          5–10m   55  (2.8×)      Fast transit → thin patches
#   Coral        10–35m   18  (0.9×)      Near-stable → cells linger
#   Grain/lagoon 35–55m   17  (0.85×)     Near-stable → open lagoon
#   Deep         55–80m   11  (0.55×)     Slow drown → background lagoon
#
# The coral zone is intentionally SLIGHTLY BELOW accommodation so cells
# aggrade slowly and spend several Myr recording coral facies before reaching
# the ooid zone. The spatial heterogeneity of the bathymetry then produces
# the patchiness: only the topographic highs (5-10m, ~7% of cells) ever
# enter the ooid zone; the rest record coral or lagoon.
#
# MODIFIER — ooid×0.2 only during 150–152.1 Ma coral bloom (t=4.8–6.9 Myr):
# The coral×5 modifier used previously caused catastrophic acceleration of
# the coral zone (69 m/Myr >> accommodation) blasting all coral-zone cells
# into the ooid zone within 0.2 Myr. Removed. The coral bloom interval is
# represented by ooid suppression (×0.2) alone, which is geologically
# appropriate — the bloom suppressed ooid factories.
# ─────────────────────────────────────────────────────────────────────────────

const MODIFIER_T_RANGE = (4.8u"Myr", 6.9u"Myr")

facies = [

   # Facies 1 — ooids
    # Narrowed to the shallowest high-energy crests and reduced in effective volume.
    # This should suppress continuous shoal blankets while preserving local shoal patches.
    # Facies 1 — ooids
# Ooids now require compact positive relief above a restricted-lagoon background.
# The window is strong but narrow: patches, not continuous belts.
WithoutCA.Facies(
    production = MultiplyProduction(
        InterpolatedProduction(
            maximum_production = 570.0u"m/Myr",
            depth_knots = [ 0.0u"m",  2.8u"m",  3.8u"m",  4.8u"m",
                            6.2u"m",  7.4u"m",  8.5u"m",  9.5u"m",
                           12.0u"m", 20.0u"m", 40.0u"m",100.0u"m",
                          600.0u"m"],
            multipliers = [  0.000,  0.000,  0.040,  0.120,
                              0.135,  0.075,  0.020,  0.000,
                              0.000,  0.000,  0.000,  0.000,
                              0.000],
        ),
        0.35; t_range = MODIFIER_T_RANGE,
    ),
    transport_coefficient = 32.0u"m/Myr",
    name = "ooids",
),

    # Facies 2 — corals
    # Shifted upward and strengthened in the 8–18 m window so corals occur more
    # commonly on the upper platform. Very low transport keeps them patchy.
    WithoutCA.Facies(
        production = InterpolatedProduction(
            maximum_production = 820.0u"m/Myr",
            depth_knots = [ 0.0u"m",  6.0u"m",  8.0u"m", 10.5u"m",
                        14.0u"m", 18.0u"m", 22.0u"m", 26.0u"m",
                        30.0u"m", 35.0u"m", 42.0u"m", 50.0u"m",
                        56.0u"m", 70.0u"m",100.0u"m",600.0u"m"],
            multipliers = [  0.000,  0.000,  0.018,  0.040,
                            0.046,  0.032,  0.012,  0.004,
                            0.010,  0.024,  0.020,  0.006,
                            0.000,  0.000,  0.000,  0.000],
        ),
        transport_coefficient = 0.8u"m/Myr",
        name = "corals",
    ),

    # Facies 3 — mud
    # Broad restricted/open lagoon matrix. Increased shallow mud support
    # makes restricted lagoon more likely where ooids/corals do not dominate.
    WithoutCA.Facies(
        production = InterpolatedProduction(
            maximum_production = 150.0u"m/Myr",
            depth_knots = [ 0.0u"m",  3.0u"m",  6.0u"m", 10.0u"m",
                           15.0u"m", 24.0u"m", 35.0u"m", 50.0u"m",
                           65.0u"m", 80.0u"m",120.0u"m",600.0u"m"],
            multipliers = [  0.16,   0.16,   0.15,   0.14,
                              0.12,   0.10,   0.08,   0.06,
                              0.04,   0.02,   0.00,   0.00],
        ),
        transport_coefficient = 18.0u"m/Myr",
        name = "mud",
    ),

    # Facies 4 — grains
    # Background peloidal/bioclastic grains, strongest in semi-open lagoon
    # and weak enough not to erase restricted lagoon.
    WithoutCA.Facies(
        production = InterpolatedProduction(
            maximum_production = 180.0u"m/Myr",
            depth_knots = [ 0.0u"m",  4.0u"m",  8.0u"m", 14.0u"m",
                           22.0u"m", 32.0u"m", 45.0u"m", 60.0u"m",
                           75.0u"m",100.0u"m",200.0u"m",600.0u"m"],
            multipliers = [  0.005,  0.012,  0.018,  0.028,
                              0.040,  0.050,  0.045,  0.030,
                              0.010,  0.000,  0.000,  0.000],
        ),
        transport_coefficient = 30.0u"m/Myr",
        name = "grains",
    ),

    # Facies 5 — sabkha / very shallow restricted cap
    # Kept minor, but inside the requested 20–100 m/Myr maximum range.
    WithoutCA.Facies(
        production = InterpolatedProduction(
            maximum_production = 25.0u"m/Myr",
            depth_knots = [ 0.0u"m",  1.0u"m",  2.0u"m",  4.0u"m",
                            6.0u"m", 10.0u"m", 20.0u"m", 50.0u"m",
                          100.0u"m",200.0u"m",600.0u"m"],
            multipliers = [  0.00,   0.08,   0.06,   0.02,
                              0.00,   0.00,   0.00,   0.00,
                              0.00,   0.00,   0.00],
        ),
        transport_coefficient = 1.0u"m/Myr",
        name = "sabkha",
    ),
]

# ─────────────────────────────────────────────────────────────────────────────
# 5. Model input
# ─────────────────────────────────────────────────────────────────────────────

input = WithoutCA.Input(
    tag                     = "eclepens-upper-malm-withoutca-patchy",
    box                     = CarboKitten.Box{Coast}(
                                  grid_size  = GRID_SIZE,
                                  phys_scale = 200.0u"m"),
    time                    = TimeProperties(
                                  t0    = 0.0u"Myr",
                                  Δt    = 0.1u"Myr",
                                  steps = 123),
    output                  = Dict(
                                  :topography => OutputSpec(
                                      slice          = (:, :),
                                      write_interval = 5),
                                  :profile    => OutputSpec(
                                      slice          = (:, div(GRID_SIZE[2], 2)),
                                      write_interval = 5)),
    initial_topography      = initial_topo,
    sea_level               = sea_level,

    # Inside 16–40 m/Myr. Slightly higher than the previous 20 m/Myr to
    # preserve shallow lagoon accommodation while allowing shoal/reef nuclei.
    subsidence_rate         = 24.0u"m/Myr",

    # Low-to-moderate reworking. Coral diffusivity is kept very low above,
    # so coral patches do not smear laterally into blankets.
    disintegration_rate     = 6.0u"m/Myr",
    lithification_time      = 2000.0u"yr",

    insolation              = 400.0u"W/m^2",
    sediment_buffer_size    = 500,
    depositional_resolution = 1.0u"m",
    facies                  = facies)

mkpath(dirname(OUTPUT_FILE))
mkpath(FIG_DIR)

@info "Running WithoutCA model..."
run_model(Model{WithoutCA}, input, OUTPUT_FILE)
@info "Done → $OUTPUT_FILE"

# ─────────────────────────────────────────────────────────────────────────────
# 6. Classification rules
# ─────────────────────────────────────────────────────────────────────────────

rules = [
    FaciesRule(
    name               = "ooid shoal",
    sediment_fractions = Dict(1 => (0.32, 1.0)),
    depth_range        = (0.0u"m", 9.0u"m"),
    wave_energy_range  = (7500.0u"W/m", Inf*u"W/m")),

    FaciesRule(
        name               = "coral patch reef",
        sediment_fractions = Dict(2 => (0.26, 1.0)),
        depth_range        = (8.0u"m", 56.0u"m"),
        wave_energy_range  = (0.0u"W/m", Inf*u"W/m")),

    FaciesRule(
        name               = "sabkha",
        sediment_fractions = Dict(5 => (0.12, 1.0)),
        depth_range        = (0.0u"m", 4.0u"m"),
        wave_energy_range  = (0.0u"W/m", 3000.0u"W/m")),

    # Expanded from 0–8 m to 0–16 m.
    # This captures the shallow protected matrix around ooid highs and between
    # coral knobs instead of pushing those cells into fallback/open lagoon.
    FaciesRule(
        name               = "restricted lagoon",
        depth_range        = (0.0u"m", 16.0u"m"),
        wave_energy_range  = (0.0u"W/m", Inf*u"W/m")),

    FaciesRule(
        name               = "semi-closed lagoon",
        depth_range        = (16.0u"m", 32.0u"m"),
        wave_energy_range  = (0.0u"W/m", Inf*u"W/m")),

    FaciesRule(
        name               = "open lagoon",
        depth_range        = (32.0u"m", 100.0u"m"),
        wave_energy_range  = (0.0u"W/m", Inf*u"W/m")),
]

# ─────────────────────────────────────────────────────────────────────────────
# 7. Read, classify
# ─────────────────────────────────────────────────────────────────────────────

header, vol  = read_volume(OUTPUT_FILE, :topography)
_,      prof = read_slice(OUTPUT_FILE,  :profile)
nx, ny       = header.grid_size

slices = vcat(
    [vol[i, :] for i in [div(nx,4), div(nx,2), 3*div(nx,4)] if i >= 1],
    [vol[:, j] for j in [div(ny,4), div(ny,2), 3*div(ny,4)] if j >= 1])

n_prof  = size(prof.sediment_thickness, 1)
sc_locs = [div(n_prof,4), div(n_prof,2), 3*div(n_prof,4)]

prod_colors = Makie.wong_colors()[1:length(facies)]
prod_labels = [f.name for f in facies]
n_cls       = length(rules) + 1
cls_colors  = Makie.wong_colors()[1:n_cls]
cls_labels  = [[r.name for r in rules]; "fallback"]

@info "Classifying..."
cls_header, cls_vol = reclassify_volume(header, vol, rules; wave_field=wave_field)
save_classified(OUTPUT_CLS_FILE, cls_header, cls_vol)

classified     = [reclassify_data(header, s, rules; wave_field=wave_field) for s in slices]
new_header     = classified[1][1]
slices_cls     = [c[2] for c in classified]
cls_prof_pairs = [reclassify_data(header, prof[loc], rules; wave_field=wave_field)
                  for loc in sc_locs]

# ─────────────────────────────────────────────────────────────────────────────
# 8. Fence diagrams
# ─────────────────────────────────────────────────────────────────────────────

function named_fence(header, slices, colors, labels; kwargs...)
    n_f = length(labels)
    fig = Figure(size=(1400, 650))
    ax  = Axis3(fig[1, 1])
    fence_diagram!(ax, header, slices;
        colormap = cgrad(colors[1:n_f], n_f; categorical=true), kwargs...)
    Legend(fig[1, 2], [PolyElement(color=colors[i]) for i in 1:n_f],
           labels; framevisible=false)
    return fig
end

fig_p = named_fence(header, slices, prod_colors, prod_labels;
    color_by=:facies, show_unconformities=10,
    show_coeval_lines=(1,5), show_sealevel=true)
save(joinpath(FIG_DIR, "fence_production.png"), fig_p; px_per_unit=2)

fig_c = named_fence(new_header, slices_cls, cls_colors, cls_labels;
    color_by=:facies, show_unconformities=10,
    show_coeval_lines=(1,5), show_sealevel=true)
save(joinpath(FIG_DIR, "fence_classified.png"), fig_c; px_per_unit=2)

# ─────────────────────────────────────────────────────────────────────────────
# 9. Map views
# ─────────────────────────────────────────────────────────────────────────────

n_frames     = length(header.axes.t[1:vol.write_interval:end])
n_cls_frames = length(cls_header.axes.t[1:cls_vol.write_interval:end])

function named_map(header, vol, colors, labels, nf; kwargs...)
    n_f = length(labels)
    fig = map_view(header, vol; times=[nf], colorbar=false,
        colormap=cgrad(colors[1:n_f], n_f; categorical=true), kwargs...)
    Legend(fig[1, 2], [PolyElement(color=colors[i]) for i in 1:n_f],
           labels; framevisible=false)
    return fig
end

fig_mp = named_map(header, vol, prod_colors, prod_labels, n_frames;
    color_by=:facies, show=:preserved, show_shoreline=false,
    mask_emerged=true, layout=:row)
save(joinpath(FIG_DIR, "mapview_production.png"), fig_mp; px_per_unit=2)

fig_mc = named_map(cls_header, cls_vol, cls_colors, cls_labels, n_cls_frames;
    color_by=:facies, show=:preserved, show_shoreline=false,
    mask_emerged=true, layout=:row)
save(joinpath(FIG_DIR, "mapview_classified.png"), fig_mc; px_per_unit=2)

# ─────────────────────────────────────────────────────────────────────────────
# 10. Stratigraphic columns
# ─────────────────────────────────────────────────────────────────────────────

phys_scale_m = ustrip(u"m", input.box.phys_scale)

fig_sp = Figure(size=(300*length(sc_locs)+120, 600))
for (k, loc) in enumerate(sc_locs)
    ax = Axis(fig_sp[1, k]; title="x = $(round(Int, loc*phys_scale_m)) m",
              xlabel="depth [m]", yreversed=false)
    stratigraphic_column!(ax, header, prof[loc]; color=prod_colors)
    k > 1 && hideydecorations!(ax; ticks=false)
end
Legend(fig_sp[1, end+1], [PolyElement(color=prod_colors[i]) for i in 1:length(facies)],
       prod_labels, "Factory"; framevisible=false)
save(joinpath(FIG_DIR, "strat_columns_production.png"), fig_sp; px_per_unit=2)

fig_sc_cls = Figure(size=(300*length(sc_locs)+150, 600))
for (k, loc) in enumerate(sc_locs)
    cls_hdr, cls_col = cls_prof_pairs[k]
    ax = Axis(fig_sc_cls[1, k]; title="x = $(round(Int, loc*phys_scale_m)) m",
              xlabel="depth [m]", yreversed=false)
    stratigraphic_column!(ax, cls_hdr, cls_col; color=cls_colors)
    k > 1 && hideydecorations!(ax; ticks=false)
end
Legend(fig_sc_cls[1, end+1], [PolyElement(color=cls_colors[i]) for i in 1:n_cls],
       cls_labels, "Environment"; framevisible=false)
save(joinpath(FIG_DIR, "strat_columns_classified.png"), fig_sc_cls; px_per_unit=2)

@info "Done. Figures in $(FIG_DIR)/"

end  # module Script
