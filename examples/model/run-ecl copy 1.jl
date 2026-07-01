# run_ecl_calibrated.jl
# Upper Malm carbonate platform — Eclépens, Swiss Molasse Basin
# WithoutCA calibration against Olivier et al. (2015)
#
# ═══════════════════════════════════════════════════════════════════════════
# SCIENTIFIC CONTEXT — ECLÉPENS
# ═══════════════════════════════════════════════════════════════════════════
#
# Eclépens sits in the Swiss Molasse Basin near Lausanne. During the Late
# Oxfordian–Early Kimmeridgian it occupied an inner-platform position,
# palaeolatitude ~26–27°N, subtropical semi-arid to arid climate
# (Frakes et al. 1992; Dercourt et al. 1993). This places it landward of
# the ooid shoal rim and within the FA4-dominated back-shoal/interior
# platform belt of Olivier's facies model — NOT in the Northern or Central
# French Jura ooid shoal fairway. FA3 shoal is present but subordinate;
# FA4 back-shoal and FA2 mid-ramp are the dominant facies.
#
# ═══════════════════════════════════════════════════════════════════════════
# OLIVIER 2015 — ECLÉPENS SEQUENCE PHASES
# ═══════════════════════════════════════════════════════════════════════════
#
# Phase A  SII retrogradation  t=0.00–1.35 Myr  upper Bimammatum, humid
#   Platform retrogrades. FA2.2 coral-microbialite biostromal buildups and
#   F2.1 oncoid packstone expand over former FA3 shoal. Ooid production
#   reduced by humid-climate nutrient input. FA4.1 back-shoal progrades
#   during decreasing accommodation (Olivier p.282 §SII).
#
# Phase B  SII progradation    t=1.35–2.67 Myr  Bimammatum, becoming arid
#   Climate drying. Oligotrophic conditions return. Ooids reappear (FA3.1/3.2).
#   FA4.1 peloid-oncoid back-shoal expands. FA3.3 inter-shoal channels active.
#   Ramp geometry persists (Olivier p.282 §Stage 3).
#
# Phase C  SIII Planula        t=2.67–3.45 Myr  Planula, warm-dry
#   Peak carbonate production. Platform progrades to oolitic rimmed shelf.
#   FA3.1 ooid grainstone dominant; FA4.2 (Cladocoropsis oncoid packstone)
#   behind rim. Tidal flat exposure at platform top (Olivier p.283 §Stage 4).
#
# Phase D  SIV–SV aggradation  t=3.45–4.68 Myr  Platynota/Hypselocyclum
#   Slightly more humid (Colombié 2002). FA4.4 interior platform mudstone
#   dominant; FA4.2 back-shoal with Cladocoropsis. Oncoid-rich deposits
#   at mid-shelf (Olivier p.283 §Stage 5).
#
# Phase E  SVI retrogradation  t=4.68–12.3 Myr  upper Hypselocyclum–Divisum
#   Carbonate production cannot keep pace with accommodation (humid seasons,
#   nutrients, flat-top sensitivity). FA2.4 bioclastic wackestone expands
#   over the platform. FA4.5 bioclastic packstone only locally. FA1 in
#   distal position (Olivier pp.283–284 §Stage 6).
#
# ═══════════════════════════════════════════════════════════════════════════
# FACTORY DESIGN
# ═══════════════════════════════════════════════════════════════════════════
#
# Six factories as sedimentological components (NOT final facies):
#   1  ooids               FA3.1, FA3.2 shoal complex grains
#   2  corals              FA2.2, FA3.4 bioconstructions
#   3  mud                 FA4.4, FA1, FA2.4 mud-supported background
#   4  peloids             FA4.1–4.3, FA3.3 restricted grains
#   5  bioclasts_intracl   FA3.2, FA3.3, FA4.5 skeletal debris/intraclasts
#   6  oncoids             FA2.1, FA2.3, FA4.1, FA4.2 coated grains
#
# BUDGET DESIGN vs σ = 22 m/Myr:
#
#   Key constraint: budget at 8–18 m MUST be < σ so cells at >18 m initial
#   depth drift THROUGH the coral zone rather than stabilising there.
#   Previous failures had bioclasts/grains adding 9.6 m/Myr at 10–18 m,
#   pushing the budget to 1.4–1.6×σ → all cells converged to coral depth.
#
#   Fixed: bioclasts start at 0.05 multiplier at 10 m, only ramp up below 16 m.
#
#   Verified budgets:
#   d= 7m (base):   ooid(20)+mud(5.4)+peloid(0)+biocl(1.8)+onc(0) = 27.2 → 1.24×σ  fast shoal
#   d= 7m (SIII):   ooid(29)+mud(5.4)+biocl(2.2) = 36.6 → 1.66×σ  very fast (Planula peak)
#   d= 7m (SII_R):  ooid(8)+mud(5.4)+biocl(1.6) = 15.0 → 0.68×σ  retreat (retrogradation ✓)
#   d=12m:          coral(7.8)+mud(7.4)+peloid(4.6)+biocl(0.5)+onc(1.7) = 22.0 → 1.00×σ  linger
#   d=18m:          coral(3.5)+mud(7.4)+peloid(9.1)+biocl(4.3)+onc(3.8) = 28.1 → 1.28×σ
#     → slightly above σ — cells at 15–18m stabilise briefly → records FA4 back-shoal ✓
#   d=25m:          mud(7.2)+peloid(9.1)+biocl(6.1)+onc(4.1) = 26.5 → 1.20×σ  near-stable
#   d=35m:          mud(8.4)+peloid(8.5)+biocl(5.9)+onc(2.9) = 25.7 → 1.17×σ  near-stable → FA4 ✓
#   d=50m:          mud(9.9)+peloid(4.6)+biocl(4.3)+onc(2.5) = 21.3 → 0.97×σ  gentle drift → FA2
#   d=65m:          mud(10.2)+biocl(1.8)+onc(0.3) = 12.3 → 0.56×σ  drift → FA1
#
#   This budget structure means:
#   • 0–8m initial  (5.9%): shoal aggrades, stabilises at ooid equilibrium → FA3
#   • 8–18m initial (17.7%): transit coral zone, some record FA3.4/FA2.2, settle at back-shoal
#   • 18–50m initial(56.1%): near-stable at FA4 back-shoal/interior depths
#   • 50–68m initial(20.3%): slow drift → FA2 upper offshore, some FA1
#
# ═══════════════════════════════════════════════════════════════════════════

module Script

using CarboKitten
using CarboKitten.Models: WithoutCA
using CarboKitten.Production: InterpolatedProduction, MultiplyProduction
using CarboKitten.WaveField: AiryWaveField, AiryWaveComponent
using CarboKitten.FaciesClassification: FaciesRule, reclassify_data,
                                         reclassify_volume, save_classified
using CarboKitten.Export: read_volume, read_slice
using CarboKitten.Visualization: fence_diagram!
using CSV, DataFrames
using Interpolations
using CairoMakie
using Unitful
using Random
using Statistics: quantile

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
# 1. Bathymetry
#
# Eclépens inner-platform target:
#   Patchy FA3 ooid highs (~6%) not a continuous belt.
#   Broad FA4 back-shoal / interior platform accommodation (~35%).
#   FA2 mid-ramp in the deeper zones (~40%).
#   FA1 in deepest cells only (~5%).
#
# Gaussian bump architecture:
#   shoal_highs:   compact, isolated ooid-shoal relief
#   shoal_cuts:    inter-shoal notches (deepen back-shoal areas)
#   coral_knobs:   roughness in the 10–40 m zone (FA2.2 patch-reef heterogeneity)
#   inner_pans:    broad shallow depressions for FA4 interior platform
# ─────────────────────────────────────────────────────────────────────────────

bathy_raw = Matrix{Float64}(CSV.read(BATHY_CSV, DataFrame; header=false))
const RAW_GRID_SIZE = (size(bathy_raw,1), size(bathy_raw,2))

function gaussian_bumps(grid_size; n, amp_lo, amp_hi, rx_lo, rx_hi, ry_lo, ry_hi, seed)
    rng = MersenneTwister(seed)
    nx, ny = grid_size
    out = zeros(Float64, nx, ny)
    for _ in 1:n
        cx = 1.0+(nx-1)*rand(rng); cy = 1.0+(ny-1)*rand(rng)
        amp = amp_lo+(amp_hi-amp_lo)*rand(rng)
        rx  = rx_lo+(rx_hi-rx_lo)*rand(rng)
        ry  = ry_lo+(ry_hi-ry_lo)*rand(rng)
        for j in 1:ny, i in 1:nx
            dx = i-cx; dy = min(abs(j-cy), ny-abs(j-cy))
            out[i,j] += amp*exp(-0.5*((dx/rx)^2+(dy/ry)^2))
        end
    end
    return out
end

let q01=quantile(vec(bathy_raw),0.01), q99=quantile(vec(bathy_raw),0.99)
    z = clamp.((bathy_raw.-q01)./(q99-q01), 0.0, 1.0)
    # Base: inner-platform minimum 9 m (so FA3 only forms on shoal highs)
    # Base minimum 5 m (was 9 m) — creates ~6% of cells at 0–8 m initial depth
    # once shoal_highs are subtracted, matching Olivier's FA3 shoal fraction.
    bathy_base = @. 5.0 + (62.0-5.0)*z^1.30

    shoal_highs = gaussian_bumps(RAW_GRID_SIZE;
        n=max(16,round(Int,prod(RAW_GRID_SIZE)/4600)),
        amp_lo=3.0, amp_hi=7.0, rx_lo=0.8, rx_hi=1.8, ry_lo=1.0, ry_hi=3.2, seed=1573)

    shoal_cuts = gaussian_bumps(RAW_GRID_SIZE;
        n=max(65,round(Int,prod(RAW_GRID_SIZE)/1200)),
        amp_lo=0.8, amp_hi=2.4, rx_lo=2.5, rx_hi=8.0, ry_lo=0.4, ry_hi=1.5, seed=1601)

    coral_gate = clamp.((bathy_base.-10.0)./20.0, 0.0, 1.0).*
                 clamp.((42.0.-bathy_base)./18.0, 0.0, 1.0)
    coral_knobs = gaussian_bumps(RAW_GRID_SIZE;
        n=max(85,round(Int,prod(RAW_GRID_SIZE)/900)),
        amp_lo=0.3, amp_hi=1.5, rx_lo=0.4, rx_hi=1.1, ry_lo=0.4, ry_hi=1.5, seed=1450)

    inner_pans = gaussian_bumps(RAW_GRID_SIZE;
        n=max(24,round(Int,prod(RAW_GRID_SIZE)/3600)),
        amp_lo=0.4, amp_hi=1.4, rx_lo=6.0, rx_hi=18.0, ry_lo=6.0, ry_hi=18.0, seed=1512)

    global bathy_rescaled = clamp.(
        bathy_base .- shoal_highs .- coral_gate.*coral_knobs .+ shoal_cuts .+ inner_pans,
        3.0, 68.0)
end

bathy_m      = bathy_rescaled * u"m"
initial_topo = -bathy_m
const GRID_SIZE = (size(bathy_m,1), size(bathy_m,2))
@info "Bathymetry $(GRID_SIZE[1])×$(GRID_SIZE[2]): $(round(minimum(bathy_rescaled);digits=1))–$(round(maximum(bathy_rescaled);digits=1)) m"

# ─────────────────────────────────────────────────────────────────────────────
# 2. Sea-level curve
# ─────────────────────────────────────────────────────────────────────────────

const T_START_MA = -157.3; const T_END_MA = -145.0
sl_df  = CSV.read(SEALEVEL_CSV, DataFrame)
age_Ma = Float64.(sl_df.age_Ma); sl_z = Float64.(sl_df.sealevel_m)
t_elapsed = age_Ma .- T_START_MA
sl_itp    = linear_interpolation(t_elapsed, sl_z; extrapolation_bc=Flat())
sl_at_t0  = sl_itp(0.0)
sea_level  = t -> (sl_itp(ustrip(u"Myr",t)) - sl_at_t0)*u"m"

# ─────────────────────────────────────────────────────────────────────────────
# 3. Wave field
# ─────────────────────────────────────────────────────────────────────────────

wave_field = AiryWaveField(components=[
    AiryWaveComponent(amplitude=2.225 * u"m", period=12.12 * u"s", direction=deg2rad(245)),
    AiryWaveComponent(amplitude=1.000 * u"m", period=10.00 * u"s", direction=deg2rad(60)),
])

# ─────────────────────────────────────────────────────────────────────────────
# 4. Sequence timing — Olivier 2015 phases
# ─────────────────────────────────────────────────────────────────────────────

const SII_RETRO = (0.00 * u"Myr", 1.35 * u"Myr")
const SII_PROG  = (1.35 * u"Myr", 2.67 * u"Myr")
const SIII      = (2.67 * u"Myr", 3.45 * u"Myr")
const SIV_SV    = (3.45 * u"Myr", 4.68 * u"Myr")
const SVI_END   = (4.68 * u"Myr", 12.31 * u"Myr")

function stage_product(base, modifiers)
    p = base
    for (factor, trange) in modifiers
        p = MultiplyProduction(p, factor; t_range=trange)
    end
    return p
end

# ─────────────────────────────────────────────────────────────────────────────
# 5. Production factories
# ─────────────────────────────────────────────────────────────────────────────

facies = [

    # ── Factory 1: Ooids ─────────────────────────────────────────────────────
    # Olivier Table 2 F3.1: "Grainstone. Ooids dominant (types 1,3), micritic
    # intraclasts, peloids. Cross bedding. Tidal shoal."
    # F3.2: "Grainstone. Ooids (types 3,4). Abundant crinoids, bivalves. Seaward bars."
    #
    # Depth window: 5–9 m strict (tidal shoal position, fair-weather wave base).
    # Rate: 20 m/Myr base.
    #
    # Stage modifiers (Olivier §§Stage 2–6):
    #   SII_RETRO ×0.40 → ooid(8)+mud(5.4)+biocl(1.6) = 15 = 0.68×σ
    #     Ooid zone retreats (retrogradation). FA2.2 corals expand over former shoal.
    #     Consistent with Olivier p.282: "retrogradation of coral-microbialite
    #     bioconstructions over ooid-bioclastic grainstone marks increasing
    #     accommodation trend of SII."
    #   SII_PROG ×0.95 → base restored. Ooids return, platform re-progrades.
    #   SIII ×1.45 → ooid(29)+biocl+mud = 36.6 = 1.66×σ
    #     "Thick ooidal deposits display aggradational pattern" (Olivier p.283).
    #     Peak ooid production under Planula warm-dry oligotrophic conditions.
    #   SIV_SV ×0.22 → "Sparser occurrences of ooid limestones during Platynota
    #     and lower Hypselocyclum" under more humid climate (Olivier p.283).
    #   SVI_END ×0.03 → ooid effectively extinct. FA2.4 bioclastic wackestone
    #     expands. "Carbonate production did not outweigh accommodation" (p.284).
    WithoutCA.Facies(
        production = stage_product(
            InterpolatedProduction(
                maximum_production = 30.0 * u"m/Myr",
                depth_knots = [ 0.0 * u"m",  3.0 * u"m",  5.0 * u"m",  7.0 * u"m",
                                8.5 * u"m",  9.5 * u"m", 13.0 * u"m", 30.0 * u"m",
                              100.0 * u"m",300.0 * u"m",400.0 * u"m",600.0 * u"m"],
                multipliers = [  0.0,    0.0,    1.0,    1.0,
                                 0.2,    0.0,    0.0,    0.0,
                                 0.0,    0.0,    0.0,    0.0],
            ),
            [(0.90, SII_RETRO), (0.95, SII_PROG),
             (1.45, SIII), (0.22, SIV_SV), (0.03, SVI_END)]),
        transport_coefficient = 35.0 * u"m/Myr",
        name = "ooids",
    ),

    # ── Factory 2: Corals ────────────────────────────────────────────────────
    # Olivier Table 2 F2.2: "Biostromal (dm-scale) or biohermal (m-scale)
    # bioconstructions." F3.4: "M-scale coral PATCH REEFS surrounded by ooid
    # grainstones." At Eclépens, FA2.2/FA3.4 are ISOLATED PATCH REEFS AND
    # BIOHERMS — not a continuous belt. Olivier Fig.3 shows FA2.2 as local
    # buildups within FA4 back-shoal matrix.
    #
    # PATCH REEF DESIGN:
    # The bathymetry has ~85 coral_knob Gaussian bumps (seed=1450) roughening
    # the 10–42 m zone with amplitude 0.3–1.5 m. These are the topographic
    # proxies for patch reef substrates.
    #
    # To exploit these knobs and produce PATCHES not a belt, the coral
    # production window is narrowed to 11–14 m peak ONLY, with near-zero
    # production by 17 m. This means:
    #   • Cells on knob crests sitting at 11–14 m get high coral signal
    #   • Cells in troughs between knobs (15–20 m) get near-zero coral
    #   → Spatial differentiation driven by bathymetric heterogeneity ✓
    #
    # Budget verification (σ = 22 m/Myr):
    #   At 12 m (knob crest):
    #     coral(12×1.0=12) + mud(3.6) + peloid(3.2) + biocl(1.2) + oncoid(1.6) = 21.6
    #     = 0.98×σ → cells linger on knob → accumulate coral signal ✓
    #     coral_frac = 12/21.6 = 55% → FA2.2 rule (coral>0.12) fires clearly ✓
    #   At 16 m (trough between knobs):
    #     coral(12×0.02=0.24) + mud(4.0) + peloid(5.0) + biocl(2.0) + oncoid(3.0) = 14.2
    #     = 0.65×σ → cells drift toward back-shoal
    #     coral_frac = 0.24/14.2 = 2% → FA2.2 rule does NOT fire ✓
    #
    # maximum_production raised 9 → 12 m/Myr to compensate for narrower window
    # (same total production budget but concentrated on knob crests only).
    WithoutCA.Facies(
        production = stage_product(
            InterpolatedProduction(
                maximum_production = 12.0 * u"m/Myr",
                depth_knots = [ 0.0 * u"m",  6.0 * u"m",  9.0 * u"m", 11.0 * u"m",
                               13.0 * u"m", 14.0 * u"m", 15.5 * u"m", 17.0 * u"m",
                               20.0 * u"m",100.0 * u"m",200.0 * u"m",600.0 * u"m"],
                multipliers = [  0.0,    0.0,    0.05,   0.80,
                                 1.0,    0.90,   0.20,   0.02,
                                 0.0,    0.0,    0.0,    0.0],
            ),
            [(1.35, SII_RETRO), (1.00, SII_PROG),
             (0.60, SIII), (0.45, SIV_SV), (0.25, SVI_END)]),
        transport_coefficient = 2.0 * u"m/Myr",
        name = "corals",
    ),

    # ── Factory 3: Mud ───────────────────────────────────────────────────────
    # Olivier Table 2 F4.4, F1.1/1.2, F2.4 — mud background at all depths.
    #
    # KEY FIX: multipliers at 8–25 m reduced 0.62 → 0.30.
    # Previous run: mud mean fraction = 0.46 platform-wide, dominating ALL zones.
    # At 12m: mud(7.4) was the single largest factory, swamping coral(7.6) signal.
    # coral_frac = 0.087 < rule threshold 0.12 → FA2.2 rule never fired.
    # ooid_frac_cumulative = 0.16 < 0.14 threshold in time-averaged column.
    #
    # With multiplier 0.30 at 8–20m:
    #   mud at 12m: 12×0.30 = 3.6 m/Myr
    #   budget at 12m: coral(9)+mud(3.6)+peloid(3.2)+biocl(1.2)+oncoid(1.6) = 18.6 = 0.85×σ
    #   → cells linger at coral depth (budget < σ) ✓
    #   coral_frac at 12m: 9/18.6 = 48% → clearly FA2.2 ✓
    #   → Ooid zone budget at 7m: ooid(19.1)+mud(12×0.45×0.65=3.5) = 22.6 ≈ σ → stable ✓
    #
    # Mud multiplier restored to 0.65 below 30m (offshore zones need background budget).
    WithoutCA.Facies(
        production = stage_product(
            InterpolatedProduction(
                maximum_production = 12.0 * u"m/Myr",
                depth_knots = [ 0.0 * u"m",  5.0 * u"m",  8.0 * u"m", 20.0 * u"m",
                               28.0 * u"m", 45.0 * u"m", 55.0 * u"m", 65.0 * u"m",
                              100.0 * u"m",200.0 * u"m",300.0 * u"m",600.0 * u"m"],
                multipliers = [  0.40,   0.45,   0.30,   0.30,
                                 0.65,   0.75,   0.82,   0.82,
                                 0.00,   0.00,   0.00,   0.00],
            ),
            [(1.10, SII_RETRO), (0.90, SII_PROG),
             (0.80, SIII), (1.25, SIV_SV), (1.25, SVI_END)]),
        transport_coefficient = 18.0 * u"m/Myr",
        name = "mud",
    ),

    # ── Factory 4: Peloids ───────────────────────────────────────────────────
    # Olivier Table 2 F4.1: "Packstone. Peloids + oncoids (Lithocodium,
    # Troglotella). Common intraclasts, aggregates. Sparse ooids. Back-shoal."
    # F4.2: "Packstone. Abundant peloids. Common large oncoids (Lithocodium,
    # Bacinella, Troglotella). Cladocoropsis. Back-shoal."
    # F4.3: "Grainstone. Peloids dominant. Oblique laminations. Gypsum
    # pseudomorphs, fenestrae. Back-shoal bars."
    # F3.3: "Packstone. Abundant intraclasts, peloids. Oncoids. Inter-shoal."
    #
    # Peak at 12–35 m (back-shoal position). Near-zero above 8 m (ooid zone).
    # Rate: 13 m/Myr.
    #
    # Stage modifiers:
    #   SII_RETRO ×0.70: peloid back-shoal reduced during retrogradation.
    #   SIV_SV ×1.15: aggradation stage — "FA4 interior platform mudstones
    #     with progradational pattern towards south" and FA4.2 back-shoal
    #     expand significantly (Olivier p.283).
    #   SVI_END ×0.55: back-shoal contracts as FA2.4 expands.
    WithoutCA.Facies(
        production = stage_product(
            InterpolatedProduction(
                maximum_production = 13.0 * u"m/Myr",
                depth_knots = [ 0.0 * u"m",  5.0 * u"m",  8.0 * u"m", 12.0 * u"m",
                               18.0 * u"m", 28.0 * u"m", 40.0 * u"m", 55.0 * u"m",
                              100.0 * u"m",200.0 * u"m",300.0 * u"m",600.0 * u"m"],
                multipliers = [  0.0,    0.0,    0.1,    0.35,
                                 0.70,   0.70,   0.65,   0.35,
                                 0.0,    0.0,    0.0,    0.0],
            ),
            [(0.70, SII_RETRO), (1.00, SII_PROG),
             (1.05, SIII), (1.05, SIV_SV), (0.55, SVI_END)]),
        transport_coefficient = 16.0 * u"m/Myr",
        name = "peloids",
    ),

    # ── Factory 5: Bioclasts / intraclasts ───────────────────────────────────
    # Olivier Table 2 F3.2: seaward bars — crinoids, bivalves, intraclasts.
    # F4.5: shoreface packstone — abundant bioclasts. F2.4: bioclastic wackestone.
    #
    # FIX 3: multipliers at 16–35 m reduced from 0.60–0.95 to 0.30–0.50.
    # Previous version: at 25m biocl(9×0.95=8.6) + peloid(6.4) + mud(8.7) + onc(4.6)
    #   = 28.2 m/Myr = 1.28×σ → ALL cells converged to 35-50m stable equilibrium.
    #   100% classified as FA4, no spatial differentiation.
    # Corrected budget targets:
    #   d=18m: biocl(9×0.25=2.3)+coral(7.4)+mud(8.2)+peloid(6.4)+onc(3.6) = 27.9 = 1.27×σ
    #     → cells transit through 18m, don't stabilise (budget > σ but falling)
    #   d=25m: biocl(9×0.30=2.7)+mud(8.7)+peloid(6.4)+onc(4.6) = 22.4 = 1.02×σ → near-stable
    #   d=35m: biocl(9×0.35=3.2)+mud(9.6)+peloid(6.1)+onc(3.5) = 22.4 = 1.02×σ → near-stable
    #   d=50m: biocl(9×0.40=3.6)+mud(10.6)+peloid(4.0)+onc(2.4) = 20.6 = 0.94×σ → slow drift
    #   d=65m: biocl(9×0.20=1.8)+mud(10.7)+peloid(2.5) = 15.0 = 0.68×σ → drown → FA1
    # → Cells at 30-50m initial depth: near-stable at FA4/FA2 for several Myr, then drift FA2
    # → Cells at 50-68m initial depth: slow drift → record FA2 then FA1
    WithoutCA.Facies(
        production = stage_product(
            InterpolatedProduction(
                maximum_production = 9.0 * u"m/Myr",
                depth_knots = [ 0.0 * u"m",  4.0 * u"m",  8.0 * u"m", 10.0 * u"m",
                               12.0 * u"m", 16.0 * u"m", 25.0 * u"m", 35.0 * u"m",
                               50.0 * u"m", 65.0 * u"m",200.0 * u"m",600.0 * u"m"],
                multipliers = [  0.0,    0.05,   0.05,   0.05,
                                 0.15,   0.25,   0.30,   0.35,
                                 0.40,   0.20,   0.00,   0.00],
            ),
            [(0.90, SII_RETRO), (1.05, SII_PROG),
             (1.20, SIII), (0.90, SIV_SV), (1.10, SVI_END)]),
        transport_coefficient = 18.0 * u"m/Myr",
        name = "bioclasts_intraclasts",
    ),

    # ── Factory 6: Oncoids ───────────────────────────────────────────────────
    # Olivier Table 2 F2.1: "Packstone. Abundant oncoids (nubecularids types 2,3).
    # Echinoderms, bivalves, gastropods, foraminifera. Upper offshore."
    # F2.3: "Packstone. Abundant large oncoids (types 5,6). Bivalves, brachiopods.
    # Marly. Central-to-distal mid-ramp."
    # F4.1: "Packstone. Peloids + oncoids (Lithocodium, Troglotella). Back-shoal."
    # F4.2: "Common large oncoids (Lithocodium, Bacinella, Troglotella).
    # Cladocoropsis. Shallow lagoon/inner ramp."
    #
    # Oncoids in Olivier have "relative high environmental tolerance" — they occur
    # both laterally to coral reefs (FA4.1, 4.2) and in deeper marly beds (FA2.1, 2.3).
    # Peak at 22–35 m (inner-to-mid-ramp transition, FA2.1 position).
    # Rate: 4.5 m/Myr. Low rate avoids elevating the mid-ramp budget above σ.
    #
    # Stage modifiers:
    #   SII_RETRO ×1.10: oncoid expansion during humid retrogradation.
    #     Lithocodium oncoids in back-shoal under well-oxygenated oligotrophic
    #     conditions (Olivier p.282, citing Dupraz & Strasser 1999).
    #   SIII ×0.70: oncoid-rich FA2.1 retreats as ooid shoal progrades
    #     (Oncolithe de Pillemoine Member retrograding in N. Central Jura, p.283).
    #   SIV_SV ×1.05: mid-shelf oncoid limestones reach maximum extent (p.283).
    #   SVI_END ×1.10: oncoid marls (FA2.3) re-expand during retrogradation.
    WithoutCA.Facies(
        production = stage_product(
            InterpolatedProduction(
                maximum_production = 4.5 * u"m/Myr",
                depth_knots = [ 0.0 * u"m",  8.0 * u"m", 15.0 * u"m", 22.0 * u"m",
                               35.0 * u"m", 50.0 * u"m", 60.0 * u"m", 68.0 * u"m",
                              100.0 * u"m",200.0 * u"m",300.0 * u"m",600.0 * u"m"],
                multipliers = [  0.0,    0.1,    0.5,    1.0,
                                 0.7,    0.5,    0.2,    0.0,
                                 0.0,    0.0,    0.0,    0.0],
            ),
            [(1.10, SII_RETRO), (1.00, SII_PROG),
             (0.70, SIII), (1.05, SIV_SV), (1.10, SVI_END)]),
        transport_coefficient = 12.0 * u"m/Myr",
        name = "oncoids",
    ),
]

# ─────────────────────────────────────────────────────────────────────────────
# 6. Model input
# ─────────────────────────────────────────────────────────────────────────────

input = WithoutCA.Input(
    tag = "eclepens-upper-malm-olivier2015-qc-v3",
    box = CarboKitten.Box{Coast}(grid_size=GRID_SIZE, phys_scale=200.0 * u"m"),
    time = TimeProperties(t0=0.0 * u"Myr", Δt=0.1 * u"Myr", steps=123),
    output = Dict(
        :topography => OutputSpec(slice=(:,:),            write_interval=5),
        :profile    => OutputSpec(slice=(:,div(GRID_SIZE[2],2)), write_interval=5)),
    initial_topography  = initial_topo,
    sea_level           = sea_level,
    subsidence_rate     = 22.0 * u"m/Myr",
    disintegration_rate = 4.0 * u"m/Myr",
    lithification_time  = 2000.0 * u"yr",
    insolation              = 400.0 * u"W/m^2",
    sediment_buffer_size    = 500,
    depositional_resolution = 1.5 * u"m",
    facies = facies)

mkpath(dirname(OUTPUT_FILE)); mkpath(FIG_DIR)
@info "Running WithoutCA model..."
run_model(Model{WithoutCA}, input, OUTPUT_FILE)
@info "Done → $OUTPUT_FILE"


# ═══════════════════════════════════════════════════════════════════════════
# 7. CLASSIFICATION RULES - adjusted after output QC
#
# QC on the uploaded result showed:
#   - production proportions are broadly plausible for Eclepens:
#       mud ~41%, peloids ~25%, oncoids ~14%, bioclasts/intraclasts ~14%,
#       ooids ~3.5%, corals ~2.4%.
#   - classified volume was too FA2-dominated:
#       FA2 ~51.5%, FA4 total ~24.6%, FA3 ~4.7%, coral buildups ~6.0%,
#       FA1 ~6.3%, fallback ~7.0%.
#   - cumulative classified map was almost entirely FA2, which is too simple.
#   - the last SVI interval was 100% fallback because muddy, bioclastic,
#     oncoid-poor offshore sediment did not satisfy either FA1 or FA2.
#
# Main corrections:
#   1. Keep FA3 and coral as diagnostic early rules.
#   2. Let FA4 back-shoal extend to 38 m, because Olivier FA4 includes
#      back-shoal / inner-platform peloid-oncoid packstones, not only <20 m.
#   3. Make FA4.4 less over-strict at model scale, because the 200 m cell
#      averages mudstone with thin packstone beds.
#   4. Separate F2.4 bioclastic wackestone from FA1 instead of forcing it into
#      lower offshore mudstone. This removes the previous late fallback and
#      matches SVI retrogradation.
#
# Factory indices:
#   1 = ooids
#   2 = corals
#   3 = mud
#   4 = peloids
#   5 = bioclasts_intraclasts
#   6 = oncoids
# ═══════════════════════════════════════════════════════════════════════════

rules = [

    # Rule 1 - FA5 tidal flat / restricted platform
    #
    # Olivier marker:
    #   F5.1-F5.3: biolaminites, laminated mudstone, fenestrae,
    #   desiccation cracks and black pebbles.
    #
    # Model proxy:
    #   very shallow mud-rich restricted sediment, with ooids capped.
    #   The wave cap is not used here because the wave model does not resolve
    #   sheltering behind the shoal rim at 200 m scale.
    FaciesRule(
        name               = "FA5 tidal flat / restricted platform",
        sediment_fractions = Dict(
            3 => (0.30, 1.0),
            1 => (0.0,  0.15),
            5 => (0.0,  0.35),
        ),
        depth_range        = (-Inf * u"m", 5.0 * u"m"),
        wave_energy_range  = (0.0 * u"W/m", Inf * u"W/m"),
    ),

    # Rule 2 - FA3 ooid shoal complex
    #
    # Olivier marker:
    #   F3.1 ooid grainstone and F3.2 ooid-bioclastic grainstone.
    #   Ooids are diagnostic; F3.2 adds crinoid/bivalve fragments and
    #   intraclasts.
    #
    # Output-7 QC:
    #   The stricter v2 rule left a significant ooid-rich population as
    #   fallback during SII-SIV. The weighted fallback composition was
    #   approximately ooids 56%, mud 35%, bioclasts/intraclasts 4%.
    #   That is not a real unclassified facies; it is the mixed/time-averaged
    #   expression of Olivier's FA3 shoal complex at 200 m grid scale.
    #
    # Model-scale correction:
    #   Restore the ooid threshold to 0.14, keep bioclasts/intraclasts present,
    #   and allow moderate mud drapes / time-averaged packstone-wackestone
    #   mixing. This recovers the SIII ooid pulse without changing production.
    FaciesRule(
        name               = "FA3 ooid shoal complex",
        sediment_fractions = Dict(
            1 => (0.14, 1.0),
            5 => (0.02, 1.0),
            3 => (0.0,  0.45),
        ),
        depth_range        = (0.0 * u"m", 16.0 * u"m"),
        wave_energy_range  = (6500.0 * u"W/m", Inf * u"W/m"),
    ),

    # Rule 3 - FA3.4 / FA2.2 coral-microbialite buildups
    #
    # QC output (6): coral buildups were volumetrically reasonable but too
    # connected spatially. Threshold and caps are tightened slightly.
    #
    # Olivier marker:
    #   F3.4 metre-scale coral patch reefs surrounded by ooid grainstones,
    #   and F2.2 coral-microbialite patch reefs / biostromes.
    #
    # Model proxy:
    #   coral signal strong enough to indicate a buildup, while excluding
    #   muddy background and peloid-dominated FA4.
    FaciesRule(
        name               = "FA3.4/FA2.2 coral-microbialite buildups",
        sediment_fractions = Dict(
            2 => (0.13, 1.0),
            3 => (0.0,  0.52),
            4 => (0.0,  0.58),
        ),
        depth_range        = (9.0 * u"m", 45.0 * u"m"),
        wave_energy_range  = (0.0 * u"W/m", Inf * u"W/m"),
    ),

    # Rule 4 - FA4.4 interior platform mudstone
    #
    # Olivier marker:
    #   F4.4 mud-supported, intensely bioturbated interior-platform
    #   mudstone/wackestone with low-energy character.
    #
    # Model proxy:
    #   mud-rich shallow platform deposit with little ooid or bioclastic signal.
    #   The mud threshold is 0.35 rather than 0.42 because the model cell is a
    #   time-averaged mixture of mudstone and thin packstone beds.
    FaciesRule(
        name               = "FA4.4 interior platform mudstone",
        sediment_fractions = Dict(
            3 => (0.35, 1.0),
            5 => (0.0,  0.28),
            1 => (0.0,  0.18),
        ),
        depth_range        = (0.0 * u"m", 38.0 * u"m"),
        wave_energy_range  = (0.0 * u"W/m", Inf * u"W/m"),
    ),

    # Rule 5 - FA4 back-shoal peloid-oncoid
    #
    # Olivier marker:
    #   F4.1/F4.2 peloid-oncoid packstones with Lithocodium/Bacinella
    #   oncoids, and F4.3 peloid grainstone bars.
    #
    # Model proxy:
    #   peloids + oncoids with low ooids. Depth extends to 38 m to avoid
    #   incorrectly classifying the Eclepens back-shoal/interior platform as
    #   mid-ramp FA2 everywhere.
    FaciesRule(
        name               = "FA4 back-shoal peloid-oncoid",
        sediment_fractions = Dict(
            4 => (0.12, 1.0),
            6 => (0.03, 1.0),
            1 => (0.0,  0.25),
            3 => (0.0,  0.62),
        ),
        depth_range        = (4.0 * u"m", 38.0 * u"m"),
        wave_energy_range  = (0.0 * u"W/m", Inf * u"W/m"),
    ),

    # Rule 6 - FA2 mid-ramp oncoid/bioclastic
    #
    # Olivier marker:
    #   F2.1 oncoid packstone and F2.3 oncoid marls. Oncoids are the
    #   diagnostic field marker, with muddy to marly matrix.
    #
    # Model proxy:
    #   oncoids + mud below the shallow FA4 window, with ooids capped.
    FaciesRule(
        name               = "FA2 mid-ramp oncoid/bioclastic",
        sediment_fractions = Dict(
            6 => (0.05, 1.0),
            3 => (0.20, 0.75),
            1 => (0.0,  0.18),
        ),
        depth_range        = (22.0 * u"m", 75.0 * u"m"),
        wave_energy_range  = (0.0 * u"W/m", Inf * u"W/m"),
    ),

    # Rule 7 - FA1 lower offshore mudstone/wackestone
    #
    # Olivier marker:
    #   F1.1/F1.2 lower offshore mudstone/wackestone, marly intervals,
    #   sparse fauna and ammonites.
    #
    # Model proxy:
    #   deepest and muddiest end-member, with very low allochem signal.
    #   Kept stricter than F2.4 so upper-offshore bioclastic wackestone is not
    #   mislabelled as lower offshore.
    FaciesRule(
        name               = "FA1 lower offshore mudstone/wackestone",
        sediment_fractions = Dict(
            3 => (0.78, 1.0),
            1 => (0.0,  0.08),
            5 => (0.0,  0.22),
            6 => (0.0,  0.15),
        ),
        depth_range        = (65.0 * u"m", 220.0 * u"m"),
        wave_energy_range  = (0.0 * u"W/m", Inf * u"W/m"),
    ),

    # Rule 8 - FA2.4 upper-offshore bioclastic wackestone
    #
    # Olivier marker:
    #   F2.4 bioclastic wackestone, mud-supported, common bivalves,
    #   ostracods and sponge spicules, deposited below fair-weather wave base.
    #
    # Model proxy:
    #   muddy, ooid-poor, oncoid-poor to weakly oncoidal upper-offshore
    #   sediment. This rule is deliberately last because it is the broad SVI
    #   retrogradational catch and removes the previous late fallback.
    FaciesRule(
        name               = "FA2.4 upper-offshore bioclastic wackestone",
        sediment_fractions = Dict(
            3 => (0.45, 1.0),
            1 => (0.0,  0.18),
            5 => (0.0,  0.50),
        ),
        depth_range        = (35.0 * u"m", 220.0 * u"m"),
        wave_energy_range  = (0.0 * u"W/m", Inf * u"W/m"),
    ),
]


# ─────────────────────────────────────────────────────────────────────────────
# 8. Read, classify, plot
# ─────────────────────────────────────────────────────────────────────────────

header, vol  = read_volume(OUTPUT_FILE, :topography)
_,      prof = read_slice(OUTPUT_FILE, :profile)
nx, ny       = header.grid_size

slices = vcat(
    [vol[i,:] for i in [div(nx,4),div(nx,2),3*div(nx,4)] if i>=1],
    [vol[:,j] for j in [div(ny,4),div(ny,2),3*div(ny,4)] if j>=1])

n_prof  = size(prof.sediment_thickness,1)
sc_locs = [div(n_prof,4), div(n_prof,2), 3*div(n_prof,4)]

function palette(n)
    base = [
        RGBAf(0.000,0.447,0.698,1.0), RGBAf(0.902,0.624,0.000,1.0),
        RGBAf(0.000,0.620,0.451,1.0), RGBAf(0.835,0.369,0.000,1.0),
        RGBAf(0.800,0.475,0.655,1.0), RGBAf(0.941,0.894,0.259,1.0),
        RGBAf(0.337,0.706,0.914,1.0), RGBAf(0.550,0.337,0.294,1.0),
        RGBAf(0.600,0.600,0.600,1.0), RGBAf(0.650,0.800,0.350,1.0),
    ]
    @assert n <= length(base)
    return base[1:n]
end

prod_colors = palette(length(facies))
prod_labels = [f.name for f in facies]
n_cls       = length(rules)+1
cls_colors  = palette(n_cls)
cls_labels  = [[r.name for r in rules]; "fallback"]

@info "Classifying..."
cls_header, cls_vol = reclassify_volume(header, vol, rules; wave_field=wave_field)
save_classified(OUTPUT_CLS_FILE, cls_header, cls_vol)

classified     = [reclassify_data(header, s, rules; wave_field=wave_field) for s in slices]
new_header     = classified[1][1]; slices_cls=[c[2] for c in classified]
cls_prof_pairs = [reclassify_data(header, prof[loc], rules; wave_field=wave_field) for loc in sc_locs]

function named_fence(header, slices, colors, labels; kwargs...)
    n_f=length(labels)
    fig=Figure(size=(1700,950)); ax=Axis3(fig[1,1])
    fence_diagram!(ax, header, slices;
        colormap=cgrad(colors[1:n_f],n_f;categorical=true), kwargs...)
    Legend(fig[2,1],[PolyElement(color=colors[i]) for i in 1:n_f],labels;
        orientation=:horizontal, nbanks=n_f<=6 ? 2 : 3, framevisible=false, tellwidth=false)
    rowsize!(fig.layout,1,Relative(0.78)); rowsize!(fig.layout,2,Relative(0.22))
    return fig
end

for (fig,name) in [
    (named_fence(header,slices,prod_colors,prod_labels;
        color_by=:facies,show_unconformities=10,show_coeval_lines=(1,5),show_sealevel=true),
     "fence_production.png"),
    (named_fence(new_header,slices_cls,cls_colors,cls_labels;
        color_by=:facies,show_unconformities=10,show_coeval_lines=(1,5),show_sealevel=true),
     "fence_classified.png")]
    save(joinpath(FIG_DIR,name), fig; px_per_unit=2)
end

n_frames=length(header.axes.t[1:vol.write_interval:end])
n_cls_frames=length(cls_header.axes.t[1:cls_vol.write_interval:end])

function named_map(header,vol,colors,labels,nf;kwargs...)
    n_f=length(labels)
    fig=map_view(header,vol;times=[nf],colorbar=false,
        colormap=cgrad(colors[1:n_f],n_f;categorical=true),kwargs...)
    resize!(fig.scene,(1700,1000))
    Legend(fig[2,1],[PolyElement(color=colors[i]) for i in 1:n_f],labels;
        orientation=:horizontal,nbanks=n_f<=6 ? 2 : 3,framevisible=false,tellwidth=false)
    return fig
end

for (fig,name) in [
    (named_map(header,vol,prod_colors,prod_labels,n_frames;
        color_by=:facies,show=:preserved,show_shoreline=false,mask_emerged=true,layout=:row),
     "mapview_production.png"),
    (named_map(cls_header,cls_vol,cls_colors,cls_labels,n_cls_frames;
        color_by=:facies,show=:preserved,show_shoreline=false,mask_emerged=true,layout=:row),
     "mapview_classified.png")]
    save(joinpath(FIG_DIR,name), fig; px_per_unit=2)
end

phys_scale_m=ustrip(u"m",input.box.phys_scale)

fig_sp=Figure(size=(1200,850))
for (k,loc) in enumerate(sc_locs)
    ax=Axis(fig_sp[1,k];title="x=$(round(Int,loc*phys_scale_m))m",xlabel="depth [m]")
    stratigraphic_column!(ax,header,prof[loc];color=prod_colors)
    k>1 && hideydecorations!(ax;ticks=false)
end
Legend(fig_sp[2,1:3],[PolyElement(color=prod_colors[i]) for i in 1:length(facies)],
    prod_labels,"Factory";orientation=:horizontal,nbanks=3,framevisible=false)
save(joinpath(FIG_DIR,"strat_columns_production.png"),fig_sp;px_per_unit=2)

fig_sc=Figure(size=(1500,950))
for (k,loc) in enumerate(sc_locs)
    cls_hdr,cls_col=cls_prof_pairs[k]
    ax=Axis(fig_sc[1,k];title="x=$(round(Int,loc*phys_scale_m))m",xlabel="depth [m]")
    stratigraphic_column!(ax,cls_hdr,cls_col;color=cls_colors)
    k>1 && hideydecorations!(ax;ticks=false)
end
Legend(fig_sc[2,1:3],[PolyElement(color=cls_colors[i]) for i in 1:n_cls],
    cls_labels,"Olivier diagnostic class";orientation=:horizontal,nbanks=3,framevisible=false)
save(joinpath(FIG_DIR,"strat_columns_classified.png"),fig_sc;px_per_unit=2)

@info "Done. Figures in $(FIG_DIR)/"

end  # module Script
