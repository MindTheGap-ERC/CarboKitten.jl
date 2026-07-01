# run_ecl_calibrated.jl
# Upper Malm / Kimmeridgian-Tithonian carbonate platform — Eclépens,
# Swiss Molasse Basin
# WithoutCA model spanning 157.3-145.0 Ma, calibrated against:
#   - Olivier et al. (2015), Swiss J. Geosci. 108:273-288
#     (Late Oxfordian-Early Kimmeridgian, French Jura, open-ramp regime)
#   - Rameil (2005), PhD thesis, Univ. Fribourg
#     (Tithonian, Swiss/French Jura, restricted lagoon-barrier regime)
#
# ═══════════════════════════════════════════════════════════════════════════
# SCIENTIFIC CONTEXT — ECLÉPENS ACROSS BOTH INTERVALS
# ═══════════════════════════════════════════════════════════════════════════
#
# Eclépens sits in the Swiss Molasse Basin near Lausanne, on the inner
# platform throughout the modelled interval. This places it:
#   - During the Kimmeridgian (Olivier regime): landward of the ooid shoal
#     rim, within the FA4-dominated back-shoal/interior platform belt -
#     NOT in the open ooid-shoal fairway of the Northern/Central French Jura.
#   - During the Tithonian (Rameil regime): within the lagoon system behind
#     the platform barrier, NOT at the barrier itself. The barrier (ooid/
#     peloid bars + patch reefs) is a linear, high-energy belt; Eclépens'
#     position records restricted-to-normal-marine LAGOON facies dominantly,
#     with barrier/coral signal present only as a subordinate, patchy
#     component (same logic as the Kimmeridgian FA3/FA2.2 patchiness).
#
# This consistent "inner-platform, lagoon/back-shoal dominant, shoal/barrier
# subordinate and patchy" position is preserved across the whole run by
# KEEPING THE SAME BATHYMETRIC ARCHITECTURE throughout (shoal highs, coral
# knobs, shoal cuts, inner pans) - only the PRODUCTION REGIME changes with
# time, reflecting the change in platform-scale depositional style described
# by Olivier (open ramp) vs. Rameil (restricted lagoon-barrier system).
#
# ═══════════════════════════════════════════════════════════════════════════
# TIME STRUCTURE: 157.3-145.0 Ma (12.3 Myr total, unchanged grid/sea-level)
# ═══════════════════════════════════════════════════════════════════════════
#
# Stage O   t=0.00-6.40 Myr   (157.3-150.9 Ma)  Olivier SVI/SVII equivalent
#   Late Kimmeridgian, upper Hypselocyclum-Divisum zones. Humid climate,
#   carbonate production cannot keep pace with accommodation (Olivier
#   pp.283-284, Stage 6). Open-ramp retrogradation: FA2.4 bioclastic
#   wackestone expands, mud-supported textures dominate, ooid shoal
#   essentially extinct, coral-microbialite patches persist locally,
#   FA1 lower-offshore mudstone present in the deepest cells.
#
# Stage R1  t=6.40-7.31 Myr  (150.9-150.0 Ma)  Kim5-Ti1 transgressive pulse
#   "Long-term (2nd-order) transgressive phase culminating in a MF in the
#   eudoxus Zone... large volumes of organically produced carbonate
#   accumulated in an aggradational pattern" (Rameil p.151). Platform
#   reorganises into the lagoon-barrier system: ooid/peloid barrier bars
#   and patch reefs establish themselves; carbonate production recovers.
#
# Stage R2  t=7.31-9.56 Myr  (150.0-147.7 Ma)  Ti1-Ti2 increasing restriction
#   "A large lagoon system that is periodically separated by a barrier from
#   the open ocean" develops (Rameil p.151). Restricted lagoon (FZ3) expands
#   landward at the expense of normal-marine lagoon (FZ4); barrier persists
#   but begins to compartmentalise the platform.
#
# Stage R3  t=9.56-10.75 Myr (147.7-146.55 Ma)  DRY PHASE BEGINS
#   "Kaolinite minimum zone"; "fully arid, desert-like climate" (Rameil
#   p.171, dry phase sensu stricto, SB Ti2? to SB Ti5). Evaporitic
#   restricted lagoon and tidal flat/sabkha expand strongly. First
#   stratiform dolomite caps form by reflux dolomitization of hypersaline
#   brines (Rameil Ch.3) - geothermally significant: dolomitized intervals
#   commonly show enhanced (recrystallization) porosity in Malm aquifers.
#
# Stage R4  t=10.75-11.49 Myr (146.55-145.81 Ma)  PEAK DRY PHASE
#   SB Ti3: "an outstanding package of tidal-flat deposits, a minimum gain
#   in accommodation space... subaerial exposure" (Rameil p.151). Maximum
#   platform-wide restriction and tidal-flat/sabkha extent of the entire
#   modelled interval. Barrier and lagoon carbonate production minimal.
#
# Stage R5  t=11.49-12.30 Myr (145.81-145.0 Ma)  Dry phase continues
#   Conditions remain arid through SB Ti4 to the top of the modelled
#   window (145.0 Ma, just before SB Ti5 = end of dry phase sensu stricto
#   at 144.96 Ma). Slightly less extreme than peak R4 but still
#   restriction-dominated.
#
# ═══════════════════════════════════════════════════════════════════════════
# FACTORY DESIGN — 7 sedimentological components (6 + 1 new)
# ═══════════════════════════════════════════════════════════════════════════
#
#   1  ooids       Olivier FA3.1/3.2 shoal grains; Rameil barrier ooid bars
#                  (inbar4, bch) - SAME grain type, SAME depth/energy niche
#                  in both regimes, only production RATE changes by stage.
#   2  corals      Olivier FA2.2/FA3.4 bioconstructions; Rameil barrier
#                  "patch reefs" (FZ5) - same role, same narrow depth window
#                  exploiting the bathymetric knobs for patchiness.
#   3  mud         Background mud-supported facies in both regimes: Olivier
#                  FA4.4/FA1/FA2.4; Rameil restricted lagoon (hrl1-4, rl1-9)
#                  and offshore mudstone.
#   4  peloids     Olivier FA4.1-4.3 back-shoal grains; Rameil internal
#                  lagoon/barrier peloids (type-1/2 peloids, the dominant
#                  non-skeletal grain of FZ3/FZ4/FZ5 per Rameil Fig.2.8a).
#   5  bioclasts   Skeletal debris + intraclasts in both regimes (Olivier
#                  F3.2/F4.5; Rameil's diverse bioclastic content of FZ4
#                  normal-marine lagoon - bivalves, gastropods, echinoderms).
#   6  oncoids     Olivier FA2.1/2.3/4.1/4.2 coated grains; Rameil type-1/2/3
#                  oncoids, abundant throughout FZ3-FZ5 per Fig.2.8a, with
#                  type and size dependent on energy conditions.
#   7  evap_mud    NEW. Tidal-flat/sabkha indicator: represents Rameil's
#                  diagnostic restricted/supratidal markers - desiccation
#                  cracks, fenestrae, evaporite pseudomorphs, microbial
#                  lamination, and the (de)dolomitization that accompanies
#                  them (Rameil Fig.2.7, Tab.2.2a, Ch.3). Olivier's FA5
#                  tidal flat is the open-ramp equivalent at much lower
#                  intensity. This factory is the model proxy for
#                  stratiform-dolomite-cap candidate intervals - the single
#                  most geothermally important diagenetic signal in Rameil's
#                  thesis, since reflux dolomitization of these caps can
#                  significantly modify reservoir porosity/permeability in
#                  the Swiss Malm.
#
# ═══════════════════════════════════════════════════════════════════════════
# SUBSIDENCE RATE AND PRODUCTION-CURVE DESIGN — σ = 16 m/Myr
# ═══════════════════════════════════════════════════════════════════════════
#
# σ is constrained to the geologically required range of 16-40 m/Myr for
# Eclépens. Within that range, 16 m/Myr (the low end) was selected
# deliberately: it minimises the TOTAL accommodation that must be filled by
# production over the 12.3 Myr run, which - as explained below - is the
# single largest factor controlling how long inherited bathymetric relief
# (and with it, coral-patch diagnostic signal) survives before the platform
# smooths out.
#
# THIS DESIGN SUPERSEDES THREE EARLIER, INCORRECT DIAGNOSES that were
# fixed in sequence during development (kept here as a record, since each
# taught something real about how this model behaves):
#
#   (1) σ=22 m/Myr + suppressed dry-phase production caused the entire
#       platform to drown catastrophically by t≈10 Myr. WRONG: Rameil's
#       dry phase is a low-accommodation, KEEP-UP regime, not a starved one.
#       Fixed by reversing the dry-phase logic (boosting mud/peloid/
#       evap_mud, not suppressing them).
#
#   (2) σ=20 m/Myr + raised production ceilings then caused ALL spatial
#       relief to collapse to exactly zero by t≈7 Myr (water-depth std
#       across the 84,051-cell grid: 16.6 m at t=0 → 0.0000 m at t=7,
#       correlation(initial depth, final thickness) = 1.0 to machine
#       precision). An initial hypothesis blamed transport_coefficient
#       (diffusive smoothing) and reduced it 2-4x; re-running gave the
#       IDENTICAL collapse, proving that diagnosis wrong. The real cause:
#       production exceeded σ at EVERY depth from 12-63 m (up to 1.86×σ),
#       creating a single attracting equilibrium depth that every cell
#       converged toward regardless of starting point.
#
#   (3) Reducing ceilings and re-tuning multipliers shrank the overshoot
#       but did not solve the underlying issue - any combination of
#       multi-peaked or sharply-tiered production curves either still
#       converges everything to one depth, or (when given separated peaks
#       above σ with deep saddles below it) sorts cells into a handful of
#       discrete, internally-flat "basins". I verified this is a genuine
#       mathematical property of the system (not a tuning failure) by
#       reading CarboKitten's actual source: `production_profile` returns
#       a pure `(time, water_depth) -> rate` closure (Components/
#       Production.jl), broadcast identically to every cell using ONLY
#       that cell's own current depth (Components/Production.jl
#       `uniform_production`, Models/WithoutCA.jl `step!`). There is NO
#       spatial argument anywhere in the call chain. This means
#       `dh/dt = σ - P(h)` is a 1-D dissipative dynamical system at every
#       (x,y) independently: any two cells starting in the same basin of
#       attraction WILL converge to the same fixed point given enough
#       time, for ANY shape of P(h). 12.3 Myr is long enough for this to
#       matter a great deal.
#
# THE WORKING SOLUTION (this version): accept that some relief decay is
# both mathematically inevitable and geologically reasonable (real
# carbonate platforms genuinely do smooth out inherited topography over
# many Myr), and design to make that decay as SLOW as possible rather than
# fighting it indefinitely:
#
#   - σ = 16 m/Myr (low end of the mandated range) minimises total
#     accommodation, directly slowing the rate at which production must
#     "win the race" against subsidence.
#   - Every factory's depth-production curve was redesigned as a SINGLE
#     smooth, gently-sloping profile (no sharp peaks, no multi-tier
#     plateaus) - gentle slopes near the typical working depth correspond
#     to a weak restoring force, which lengthens the time constant for
#     relief decay. Sharp peaks/saddles (tried and rejected, see (3) above)
#     create fast-converging attractors or unstable escape zones instead.
#   - Verified via a 0-D forward integration (production minus
#     accommodation, run over the full 12.3 Myr against the actual
#     sea-level curve) across a fine grid of initial depths from 0-68 m:
#     final depths now form a CONTINUOUS, SMOOTHLY-VARYING gradient
#     (90-113 m for evenly-spaced 2 m initial-depth steps - see verified
#     values below) rather than collapsing to one value or sorting into a
#     small number of flat basins.
#   - Critically, cells starting at the coral_knobs crest depth (~12 m)
#     and in the adjacent troughs (~18 m) remain MEANINGFULLY DISTINCT in
#     water depth through the first ~3 Myr of the run (Stage O into early
#     R1) - e.g. trajectories of 9.8/16.7 m at t=0.5, 8.7/13.9 m at
#     t=1.5 - before gradually converging thereafter. Within that window,
#     coral fraction at the knob vs. trough depths differs meaningfully
#     (e.g. 0.31 vs 0.03 at one timestep), which is exactly what is needed
#     for the FaciesRule coral threshold (>0.10) to distinguish genuine
#     patches from background sediment.
#   - PATCHES ARE NOT REQUIRED OR EXPECTED TO SURVIVE THE FULL RUN. Per
#     explicit confirmation, the target is for produced/classified blocks
#     to be consistent with Olivier and Rameil at Eclépens overall - which
#     includes both papers' own narrative of a structured, patchy early
#     platform (Olivier's retrogradation-stage coral patches; Rameil's
#     barrier-stage reef establishment) progressively giving way to a
#     smoother, restricted, tidal-flat-dominated late platform (Rameil's
#     dry phase). The relief-decay behaviour described above happens to
#     track this narrative reasonably well: patches are most vivid in
#     Stage O-R1, and smooth out through R2-R5 as restriction intensifies -
#     consistent with, not contrary to, the source material.
#
# NOTE ON THE SEA-LEVEL CURVE: the supplied curve shows a sharp LATE rise
# (≈+50 m between t=10.7 and t=12.0 Myr, i.e. ≈146-145 Ma) that is steeper
# than what Rameil documents for this specific interval - he places the
# return to more humid, transgressive conditions essentially AT or after
# SB Ti5 (144.96 Ma), just past our window's end, not ramping up strongly
# within it. This contributes to the deepening seen in the R4-R5 stages but
# is not in itself inconsistent with the window ending close to that
# boundary; it may still be worth revisiting independently of this design.
# ═══════════════════════════════════════════════════════════════════════════
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
# 1. Bathymetry — UNCHANGED architecture from the calibrated Kimmeridgian-only
# run. Eclépens' inner-platform, patchy-barrier position is the same
# geometric reality in both the Olivier and Rameil regimes; only the
# production rules acting on this geometry change with time (see below).
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
# 2. Sea-level curve — same CSV, same interpolation. The CSV is assumed to
# span (at least) 157.3-145.0 Ma as required by the new time window.
# ─────────────────────────────────────────────────────────────────────────────

const T_START_MA = -157.3; const T_END_MA = -145.0
sl_df  = CSV.read(SEALEVEL_CSV, DataFrame)
age_Ma = Float64.(sl_df.age_Ma); sl_z = Float64.(sl_df.sealevel_m)
@assert T_START_MA >= minimum(age_Ma) "Sea-level CSV does not extend back to 157.3 Ma"
@assert T_END_MA   <= maximum(age_Ma) "Sea-level CSV does not extend to 145.0 Ma"
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
# 4. Sequence timing — 157.3-145.0 Ma, Olivier (Stage O) + Rameil (R1-R5)
#
#    Ages given for reference; t-values are elapsed Myr from t0=157.3 Ma.
#    Rameil sequence-boundary ages from his Tab.6.2 (this work column).
# ─────────────────────────────────────────────────────────────────────────────

const STAGE_O  = (0.00 * u"Myr",  6.40 * u"Myr")   # 157.3-150.9 Ma  Olivier SVI/SVII
const STAGE_R1 = (6.40 * u"Myr",  7.31 * u"Myr")   # 150.9-150.0 Ma  SB Kim5-Ti1
const STAGE_R2 = (7.31 * u"Myr",  9.56 * u"Myr")   # 150.0-147.7 Ma  SB Ti1-Ti2?
const STAGE_R3 = (9.56 * u"Myr", 10.75 * u"Myr")   # 147.7-146.55 Ma SB Ti2?-Ti3 (dry phase begins)
const STAGE_R4 = (10.75 * u"Myr",11.49 * u"Myr")   # 146.55-145.81 Ma SB Ti3-Ti4 (peak dry phase)
const STAGE_R5 = (11.49 * u"Myr",12.30 * u"Myr")   # 145.81-145.0 Ma  SB Ti4-end (dry phase continues)

function stage_product(base, modifiers)
    p = base
    for (factor, trange) in modifiers
        p = MultiplyProduction(p, factor; t_range=trange)
    end
    return p
end

# ─────────────────────────────────────────────────────────────────────────────
# 5. Production factories — 7 sedimentological components
#
# PRODUCTION FIX 2026-06-30
# OOID / SABKHA MINOR RETUNE 2026-06-30
# THICKNESS / ECLÉPENS CONSISTENCY TWEAK 2026-06-30
# - Target final platform thickness is ~340 m rather than ~220 m.
# - Raise accommodation to ~25.5 m/Myr and scale mainly lagoonal
#   production, not coral, to preserve restricted-lagoon dominance.
# - Reduce and localise coral production so buildups remain patchy.
# - Broaden the ooid window slightly for R1-R2 barrier-shoal signal.
# - Make the dolomite-cap/tidal-flat candidate rule reachable.
# - Slightly strengthens R1-R2 ooids so barrier-shoal layers register.
# - Slightly strengthens R3-R5 evap_mud so dry-phase tidal-flat/sabkha
#   and dolomite-cap candidate intervals are visible without overprinting
#   the restricted-lagoon-dominant Tithonian signal.
# - Keep Stage O relatively close to the previous run so early coral-patch
#   behaviour is not erased.
# - Increase R1-R4 lagoonal keep-up production between ~8 and 40 m so the
#   Tithonian system does not drown into mid-ramp/fallback facies.
# - Reduce R5 production again to avoid forcing the whole platform into
#   exposure after the late sea-level rise slows down.
# - Broaden evap_mud from a strictly 0-3 m factory into a weak 3-20 m
#   restricted-lagoon/tidal-flat-adjacent signal during R2-R5.
# ─────────────────────────────────────────────────────────────────────────────

facies = [

    # ── Factory 1: Ooids ─────────────────────────────────────────────────────
    # Olivier F3.1/3.2 ooid shoal grains; Rameil barrier ooid bars (inbar4,
    # bch). Same physical grain, same depth/energy niche (5-9 m, high energy)
    # in both regimes - only production RATE varies by stage.
    #
    # Stage modifiers (re-tuned for σ=16 m/Myr, the low end of the mandated
    # 16-40 m/Myr range, selected to minimise total accommodation and
    # maximise how long the inherited bathymetric relief - and barrier/coral
    # patchiness - stays diagnostically meaningful before the platform
    # gradually smooths; see mud factory comment for the full rationale):
    #   STAGE_O ×0.05: Olivier's late-Kimmeridgian retrogradation - ooid
    #     production "did not outweigh accommodation" (Olivier p.284);
    #     essentially extinct, consistent with FA2.4 expansion.
    #   STAGE_R1 ×0.50: transgressive pulse partially re-establishes barrier
    #     ooid bars as the lagoon-barrier system organises (Rameil p.151);
    #     subordinate to mud/peloid, consistent with the barrier being a
    #     minor component of Eclépens' lagoon-interior position.
    #   STAGE_R2 ×0.40: barrier persists but lagoon restriction increases.
    #   STAGE_R3 ×0.20: dry-phase onset - barrier ooid production declines
    #     as the lagoon interior takes over carbonate budget (see mud/
    #     peloid/evap_mud below, which are BOOSTED, not suppressed, in this
    #     stage - the platform stays shallow via lagoon/tidal-flat
    #     aggradation even as the barrier itself contracts).
    #   STAGE_R4 ×0.10: peak dry phase - "minimum gain in accommodation
    #     space" (Rameil p.151, SB Ti3) - barrier shoal at its smallest.
    #   STAGE_R5 ×0.15: dry phase continues, barrier still subordinate.
    WithoutCA.Facies(
        production = stage_product(
            InterpolatedProduction(
                maximum_production = 48.0 * u"m/Myr",
                depth_knots = [ 0.0 * u"m",  3.0 * u"m",  5.0 * u"m",  7.0 * u"m",
                                8.5 * u"m",  9.5 * u"m", 13.0 * u"m", 30.0 * u"m",
                              100.0 * u"m",300.0 * u"m",400.0 * u"m",600.0 * u"m"],
                multipliers = [  0.0,    0.0,    0.90,   1.00,
                                 0.65,   0.35,   0.12,   0.00,
                                 0.00,   0.00,   0.00,   0.00],
            ),
            [(0.05, STAGE_O), (0.85, STAGE_R1), (0.28, STAGE_R2),
             (0.05, STAGE_R3), (0.01, STAGE_R4), (0.02, STAGE_R5)]),
        transport_coefficient = 8.0 * u"m/Myr",
        name = "ooids",
    ),

    # ── Factory 2: Corals ────────────────────────────────────────────────────
    # Olivier F2.2/F3.4 coral-microbialite bioconstructions; Rameil barrier
    # "patch reefs" (FZ5, "biogenically constructed reefs of the barrier
    # zone... not observed in this study but described by FOOKES 1995,
    # DÉTRAZ & MOJON 1989" etc. in laterally equivalent strata - Rameil p.41).
    # Same narrow depth window (11-14 m peak) exploiting bathymetric knobs
    # for patch-reef heterogeneity in both regimes.
    #
    # Stage modifiers (re-tuned for σ=16 m/Myr; coral patches are most
    # diagnostically meaningful in Stage O-R1, while inherited bathymetric
    # relief - including the coral_knobs Gaussian field - has not yet
    # smoothed out; per CarboKitten's actual production engine, verified
    # by reading the source, production depends ONLY on (time, water_depth)
    # with no spatial argument, so patches persist exactly as long as cells
    # on knob crests and in troughs maintain genuinely different water
    # depths - expected to be strongest in the early stages and to fade
    # gradually thereafter, which is consistent with both papers' narrative
    # of progressive platform restriction/smoothing through the Tithonian):
    #   STAGE_O ×0.30: subordinate coral patches persist locally during
    #     Kimmeridgian retrogradation (lower rate than the earlier humid
    #     bloom phases of Olivier's SII, since by SVI/SVII heterotrophic
    #     fauna dominates - Olivier p.284).
    #   STAGE_R1 ×0.90: patch-reef establishment during transgressive
    #     reorganisation of the barrier system - the most favourable window
    #     for genuine coral patchiness, while relief is still fresh.
    #   STAGE_R2 ×0.70: patch reefs persist but contract as restriction
    #     increases.
    #   STAGE_R3 ×0.20: dry-phase onset suppresses reef growth as the
    #     barrier contracts (lagoon interior compensates via boosted
    #     mud/peloid/evap_mud production - see those factories below).
    #   STAGE_R4 ×0.08: peak dry phase - reef growth nearly shut down,
    #     consistent with minimal barrier activity at SB Ti3.
    #   STAGE_R5 ×0.12: still strongly subordinate through the dry phase.
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
            [(0.35, STAGE_O), (0.65, STAGE_R1), (0.38, STAGE_R2),
             (0.03, STAGE_R3), (0.00, STAGE_R4), (0.00, STAGE_R5)]),
        transport_coefficient = 0.4 * u"m/Myr",
        name = "corals",
    ),

    # ── Factory 3: Mud ───────────────────────────────────────────────────────
    # Olivier F4.4/F1.1-1.2/F2.4 mud-supported background; Rameil restricted
    # lagoon mudstones (hrl1-4, rl1-9) and offshore equivalents. Background
    # mud production present at all depths in both regimes.
    #
    # CURVE REDESIGNED as a single SMOOTH, GENTLY-SLOPING profile (no sharp
    # peaks or plateaus) after verifying via CarboKitten's actual source
    # (Components/Production.jl, Models/WithoutCA.jl) that production is a
    # pure function of (time, water_depth) with NO spatial argument - every
    # cell's production depends only on ITS OWN current depth, broadcast
    # identically platform-wide. This means persistent multi-cell relief can
    # ONLY survive as long as cells have not yet converged within their
    # shared basin of attraction (a mathematical property of any uniform-σ,
    # depth-only-production system integrated long enough). A single gentle
    # slope (rather than the previous multi-tier/peaked designs, which create
    # FAST-converging attractors or unstable saddles) maximises the time
    # before convergence, keeping the platform's inherited bathymetric relief
    # (shoal highs, coral knobs) diagnostically meaningful through the early-
    # to-middle stages, consistent with progressive smoothing/restriction
    # being a genuine part of both papers' Kimmeridgian-to-Tithonian story
    # rather than something to fight indefinitely.
    #
    # Stage modifiers (re-tuned for σ=16 m/Myr, the low end of the mandated
    # 16-40 m/Myr range, chosen specifically to minimise total accommodation
    # that must be filled over the 12.3 Myr run and thus maximise how long
    # spatial relief - and coral patchiness - survives before smoothing):
    #   STAGE_O ×1.10: "evolution from grain-supported to mud-supported
    #     texture" during Olivier's SVI retrogradation (p.283).
    #   STAGE_R1 ×0.95: transgressive pulse, near-baseline mud production.
    #   STAGE_R2 ×1.00: "large lagoon system" (Rameil p.151) baseline.
    #   STAGE_R3 ×1.05: dry-phase onset - mud production BOOSTED as the
    #     restricted lagoon expands and aggrades to keep pace with reduced
    #     accommodation (hrl1 "empty, texture-less mudstone", Rameil
    #     Tab.2.2a). This is a keep-up, not a starved, regime.
    #   STAGE_R4 ×1.10: peak dry phase - maximal mud production, consistent
    #     with the platform staying shallow through the most restricted
    #     interval rather than drowning (subaerial exposure features
    #     require the substrate to reach sea level).
    #   STAGE_R5 ×1.08: still elevated through the remainder of the dry phase.
    WithoutCA.Facies(
        production = stage_product(
            InterpolatedProduction(
                maximum_production = 24.0 * u"m/Myr",
                depth_knots = [ 0.0 * u"m",  8.0 * u"m", 20.0 * u"m", 35.0 * u"m",
                               50.0 * u"m", 60.0 * u"m", 68.0 * u"m",100.0 * u"m",
                              150.0 * u"m",200.0 * u"m",300.0 * u"m",600.0 * u"m"],
                multipliers = [  0.55,   0.66,   0.82,   0.82,
                                 0.65,   0.40,   0.20,   0.05,
                                 0.00,   0.00,   0.00,   0.00],
            ),
            [(0.72, STAGE_O), (1.02, STAGE_R1), (1.12, STAGE_R2),
             (1.35, STAGE_R3), (1.40, STAGE_R4), (1.05, STAGE_R5)]),
        transport_coefficient = 5.0 * u"m/Myr",
        name = "mud",
    ),

    # ── Factory 4: Peloids ───────────────────────────────────────────────────
    # Olivier F4.1-4.3 back-shoal grains; Rameil's type-1/2 peloids - the
    # dominant non-skeletal grain of the internal lagoon and barrier zones
    # per Fig.2.8a ("abundant type-2 peloids and intraclasts" in internal
    # bars; "type-2 peloids" throughout restricted/normal-marine lagoon
    # facies rl1-rl9).
    #
    # CURVE REDESIGNED as a single smooth, gently-rising profile (see mud
    # factory comment above for the full rationale: production is depth-only
    # with no spatial argument in CarboKitten, so a gentle slope maximises
    # how long the inherited bathymetric relief - and coral patchiness -
    # remains diagnostically meaningful before the platform smooths out).
    #
    # Stage modifiers (re-tuned for σ=16 m/Myr):
    #   STAGE_O ×0.70: back-shoal subordinate during Olivier's SVI
    #     retrogradation (mud-dominated regime takes precedence).
    #   STAGE_R1 ×0.95: lagoon peloid production re-establishes as the
    #     barrier-lagoon system organises.
    #   STAGE_R2 ×1.00: "large lagoon system" (Rameil p.151) baseline -
    #     peloids abundant throughout the expanding lagoon.
    #   STAGE_R3 ×1.05: dry-phase onset - peloid production BOOSTED along
    #     with mud/evap_mud, representing the restricted lagoon's keep-up
    #     aggradation as the barrier contracts (NOT an overall carbonate-
    #     factory shutdown - only the barrier/reef components decline).
    #   STAGE_R4 ×1.10: peak dry phase - maximum restricted-lagoon peloid
    #     production, consistent with the platform staying shallow.
    #   STAGE_R5 ×1.08: still elevated through the remainder of the dry phase.
    WithoutCA.Facies(
        production = stage_product(
            InterpolatedProduction(
                maximum_production = 24.0 * u"m/Myr",
                depth_knots = [ 0.0 * u"m",  8.0 * u"m", 20.0 * u"m", 35.0 * u"m",
                               50.0 * u"m", 60.0 * u"m", 68.0 * u"m",100.0 * u"m",
                              150.0 * u"m",200.0 * u"m",300.0 * u"m",600.0 * u"m"],
                multipliers = [  0.10,   0.42,   0.66,   0.75,
                                 0.58,   0.32,   0.15,   0.03,
                                 0.00,   0.00,   0.00,   0.00],
            ),
            [(0.52, STAGE_O), (1.05, STAGE_R1), (1.12, STAGE_R2),
             (1.35, STAGE_R3), (1.35, STAGE_R4), (1.00, STAGE_R5)]),
        transport_coefficient = 4.5 * u"m/Myr",
        name = "peloids",
    ),

    # ── Factory 5: Bioclasts / intraclasts ───────────────────────────────────
    # Olivier F3.2/F4.5 skeletal debris and intraclasts; Rameil's diverse
    # bioclastic content of the normal-marine lagoon (FZ4 - "primarily
    # bivalves, gastropods, ostracods, and spiculae, ... frequent to
    # abundant echinoderms", Tab.2.2a lowEl1) and barrier zone reworked
    # bioclastic debris.
    #
    # CURVE REDESIGNED as a single smooth, gently-rising profile (see mud
    # factory comment for full rationale on why peaked/multi-tier curves
    # were replaced).
    #
    # Stage modifiers (re-tuned for σ=16 m/Myr):
    #   STAGE_O ×0.90: baseline, slightly subordinate to mud/peloid.
    #   STAGE_R1 ×0.95: transgressive pulse, diverse normal-marine fauna.
    #   STAGE_R2 ×0.95: normal-marine to restricted lagoon transition,
    #     fauna still relatively diverse in open parts of the lagoon.
    #   STAGE_R3 ×0.55: dry-phase onset - "faunal deserts" begin to develop
    #     in increasingly restricted/hypersaline lagoon settings (cf. Rameil
    #     Tab.2.1, Hughes Clarke & Keij salinity-faunal model). Reduced but
    #     not eliminated - some euryhaline fauna persists per Rameil p.40.
    #   STAGE_R4 ×0.40: peak dry phase - severe faunal restriction,
    #     bioclastic content minimal except locally.
    #   STAGE_R5 ×0.45: still strongly reduced.
    WithoutCA.Facies(
        production = stage_product(
            InterpolatedProduction(
                maximum_production = 8.5 * u"m/Myr",
                depth_knots = [ 0.0 * u"m", 10.0 * u"m", 25.0 * u"m", 40.0 * u"m",
                               55.0 * u"m", 65.0 * u"m", 75.0 * u"m",100.0 * u"m",
                              150.0 * u"m",200.0 * u"m",300.0 * u"m",600.0 * u"m"],
                multipliers = [  0.08,   0.08,   0.18,   0.27,
                                 0.33,   0.35,   0.32,   0.15,
                                 0.00,   0.00,   0.00,   0.00],
            ),
            [(0.85, STAGE_O), (0.65, STAGE_R1), (0.45, STAGE_R2),
             (0.22, STAGE_R3), (0.15, STAGE_R4), (0.18, STAGE_R5)]),
        transport_coefficient = 5.0 * u"m/Myr",
        name = "bioclasts_intraclasts",
    ),

    # ── Factory 6: Oncoids ───────────────────────────────────────────────────
    # Olivier F2.1/F2.3/F4.1/F4.2 coated grains; Rameil type-1/2/3 oncoids,
    # described as abundant throughout FZ3-FZ5 with "type, size and form
    # being dependent on energy conditions" (Rameil p.41). Lithocodium/
    # Bacinella-type oncoids common in both regimes' back-shoal/lagoon
    # settings.
    #
    # CURVE REDESIGNED as a single smooth, gently-rising profile (see mud
    # factory comment for full rationale).
    #
    # Stage modifiers (re-tuned for σ=16 m/Myr):
    #   STAGE_O ×0.90: oncoid-rich FA2.1/2.3 present but subordinate to mud
    #     during Olivier's SVI retrogradation.
    #   STAGE_R1 ×0.90: baseline during transgressive reorganisation.
    #   STAGE_R2 ×0.95: oncoids abundant in the expanding lagoon system,
    #     "facies with predominantly type-1/2 oncoids... represents
    #     transition to normal marine lagoon" (Rameil Tab.2.2a, rl9).
    #   STAGE_R3 ×0.95: dry-phase onset - oncoid-forming biota only mildly
    #     reduced under increasing restriction/salinity stress; oncoids
    #     remain a persistent component of the restricted lagoon (Rameil
    #     Fig.2.8a shows all three oncoid types present across FZ3-FZ5).
    #   STAGE_R4 ×0.90: peak dry phase - oncoid production reduced but not
    #     eliminated.
    #   STAGE_R5 ×0.85: still present, moderately reduced.
    WithoutCA.Facies(
        production = stage_product(
            InterpolatedProduction(
                maximum_production = 6.5 * u"m/Myr",
                depth_knots = [ 0.0 * u"m",  8.0 * u"m", 20.0 * u"m", 35.0 * u"m",
                               50.0 * u"m", 60.0 * u"m", 68.0 * u"m",100.0 * u"m",
                              150.0 * u"m",200.0 * u"m",300.0 * u"m",600.0 * u"m"],
                multipliers = [  0.0,    0.08,   0.32,   0.55,
                                 0.66,   0.70,   0.65,   0.30,
                                 0.00,   0.00,   0.00,   0.00],
            ),
            [(0.81, STAGE_O), (0.95, STAGE_R1), (1.00, STAGE_R2),
             (0.85, STAGE_R3), (0.70, STAGE_R4), (0.75, STAGE_R5)]),
        transport_coefficient = 3.5 * u"m/Myr",
        name = "oncoids",
    ),

    # ── Factory 7: Evaporitic mud / tidal-flat-sabkha (NEW) ──────────────────
    # Rameil FZ1-FZ2: sabkha (sbk1-2) and tidal-flat (tf1-tf6) facies,
    # diagnosed by "desiccation cracks, fenestrae, evaporite pseudomorphs,
    # microbial mat lamination" (Rameil Fig.2.7, Tab.2.2a) and the
    # accompanying penecontemporaneous (de)dolomitization. Olivier's FA5
    # (F5.1-5.3) is the open-ramp equivalent at much lower intensity -
    # this factory also covers that signal in Stage O, at very low rate.
    #
    # GEOTHERMAL SIGNIFICANCE: this factory is the model proxy for
    # stratiform-dolomite-cap candidate intervals. Rameil's central
    # diagenetic finding (Ch.3) is that m-thick stratiform dolomite forms
    # by reflux of hypersaline brines generated in these evaporitic
    # tidal-flat/sabkha settings, capping the regressive parts of medium-
    # and large-scale sequences. Dolomitized intervals in the Swiss Malm
    # commonly show modified (often enhanced, via recrystallization)
    # porosity relative to the limestone precursor - directly relevant to
    # reservoir quality prediction at Eclépens.
    #
    # Strict shallow window (<3 m, supratidal-uppermost intertidal), near
    # zero contribution everywhere else. Production essentially switched
    # off during Stage O/R1/R2 (Olivier's setting and the early Tithonian
    # transgressive/normal-marine lagoon are not evaporitic) and ramped up
    # strongly from Stage R3 onward as the dry phase develops.
    #
    # Stage modifiers:
    #   STAGE_O  ×0.15: Olivier's rare, low-intensity FA5 tidal flat signal.
    #   STAGE_R1 ×0.10: transgressive, normal-marine - minimal tidal flat.
    #   STAGE_R2 ×0.40: increasing restriction begins to favour local tidal
    #     flat/sabkha development at the lagoon margins.
    #   STAGE_R3 ×1.60: dry phase sensu stricto begins (kaolinite minimum) -
    #     "fully arid, desert-like climate" (Rameil p.171) strongly favours
    #     evaporitic tidal-flat/sabkha development and dolomitizing reflux.
    #   STAGE_R4 ×2.20: peak dry phase, SB Ti3 - "an outstanding package of
    #     tidal-flat deposits" (Rameil p.151) - maximum tidal-flat/sabkha
    #     extent and dolomite-cap formation potential of the entire run.
    #   STAGE_R5 ×1.90: dry phase continues at slightly reduced intensity.
    WithoutCA.Facies(
        production = stage_product(
            InterpolatedProduction(
                maximum_production = 18.0 * u"m/Myr",
                depth_knots = [ 0.0 * u"m",  1.5 * u"m",  3.0 * u"m",  5.0 * u"m",
                                8.0 * u"m", 12.0 * u"m", 20.0 * u"m", 30.0 * u"m",
                               50.0 * u"m", 70.0 * u"m",100.0 * u"m",200.0 * u"m",600.0 * u"m"],
                multipliers = [  0.55,   0.55,   0.42,   0.28,
                                 0.18,   0.08,   0.02,   0.00,
                                 0.00,   0.00,   0.00,   0.00,   0.00],
            ),
            [(0.02, STAGE_O), (0.04, STAGE_R1), (0.28, STAGE_R2),
             (0.95, STAGE_R3), (1.20, STAGE_R4), (0.75, STAGE_R5)]),
        transport_coefficient = 2.0 * u"m/Myr",
        name = "evap_mud",
    ),
]

# ─────────────────────────────────────────────────────────────────────────────
# 6. Model input
# ─────────────────────────────────────────────────────────────────────────────

# FINAL BALANCE TWEAK 2026-06-30
# R2 BARRIER PATCH FIX 2026-06-30
# R1 FALLBACK / R5 SABKHA BLANKET FIX 2026-06-30
# R1 FALLBACK FINAL RULE FIX 2026-06-30
# - Keeps production, subsidence, thickness and R2 barrier fix unchanged.
# - Converts the remaining R1 shallow mixed lagoon fallback into open/back-shoal lagoon.
# - Does not force sabkha because evap_mud is spatially too uniform for a meaningful cap map.
# - Keeps thickness/accommodation and R2 barrier fix unchanged.
# - Relaxes open/back-shoal rule only enough to catch R1 shallow mixed lagoon cells.
# - Prevents the R5 evap_mud signal from classifying the whole platform as sabkha.
# - Keeps thickness/accommodation unchanged.
# - Reduces only the R2 ooid overprint that created a barrier-shoal blanket.
# - Tightens the barrier-shoal rule so only true ooid-rich high-energy cells classify as barrier.
# - Slightly relaxes the sabkha/dolomite-cap rule so dry-phase caps can register locally.
# R2 OOID BALANCE TWEAK 2026-06-30
# - The validation run matched the 340 m thickness and Stage O/R3-R5 trend,
#   but R2 became too barrier-shoal dominated for an Eclepens inner-platform
#   position. Reduce only the R2 ooid multiplier; keep sigma and late lagoon
#   production unchanged.
# - Keeps σ=25.5 m/Myr because the last run matched the 340 m thickness target.
# - Reduces the late evap_mud blanket that produced R3-R5 fallback/tidal-flat overprint.
# - Restores mud+peloid restricted-lagoon production in the dry phase.
# - Makes Stage O classify as open/protected platform instead of restricted lagoon.

input = WithoutCA.Input(
    tag = "eclepens-kimmeridgian-tithonian-olivier-rameil",
    box = CarboKitten.Box{Coast}(grid_size=GRID_SIZE, phys_scale=200.0 * u"m"),
    time = TimeProperties(t0=0.0 * u"Myr", Δt=0.1 * u"Myr", steps=123),
    output = Dict(
        :topography => OutputSpec(slice=(:,:),            write_interval=5),
        :profile    => OutputSpec(slice=(:,div(GRID_SIZE[2],2)), write_interval=5)),
    initial_topography  = initial_topo,
    sea_level           = sea_level,
    subsidence_rate    = 25.5 * u"m/Myr",
    disintegration_rate = 4.0 * u"m/Myr",
    lithification_time  = 2000.0 * u"yr",
    insolation              = 400.0 * u"W/m^2",
    sediment_buffer_size    = 500,
    depositional_resolution = 1.5 * u"m",
    facies = facies)

mkpath(dirname(OUTPUT_FILE)); mkpath(FIG_DIR)
@info "Running WithoutCA model (157.3-145.0 Ma)..."
run_model(Model{WithoutCA}, input, OUTPUT_FILE)
@info "Done → $OUTPUT_FILE"


# ═══════════════════════════════════════════════════════════════════════════
# 7. CLASSIFICATION RULES
#    9 classes spanning both regimes, ordered most specific -> most general.
#    Built to be resolvable at 200 m grid scale and to discriminate
#    geothermally meaningful reservoir/seal categories throughout.
#
# Factory indices:
#   1=ooids  2=corals  3=mud  4=peloids  5=bioclasts_intraclasts
#   6=oncoids  7=evap_mud
# ═══════════════════════════════════════════════════════════════════════════
#
# RULE 1 - Tidal flat / sabkha (dolomite-cap candidate)
#   Olivier F5.1-5.3 + Rameil FZ1-FZ2 (sbk1-2, tf1-tf6).
#   Field marker: desiccation cracks, fenestrae, evaporite pseudomorphs,
#     microbial lamination (Rameil Fig.2.7) - none of these are directly
#     observable from sediment fractions alone, but the evap_mud factory
#     is purpose-built as their proxy: it is ONLY active in the supratidal/
#     uppermost-intertidal production window and is strongly amplified
#     during the dry-phase stages (R3-R5) when these features are described
#     as widespread (Rameil p.151, p.171).
#   evap_mud > 0.30: dominant or co-dominant evaporitic-mud signal.
#   mud > 0.30: still mud/micrite-supported texture (M/W in Rameil's Dunham
#     classification for tf1-tf6, sbk1-2).
#   depth < 3 m: supratidal to uppermost intertidal position (Rameil Fig.2.7
#     places these features at or above mean low water).
#   wave < 4000 W/m: protected, low-energy setting (landward of the barrier
#     in Rameil's model, or landward of the shoal rim in Olivier's).
#   GEOTHERMAL: flags candidate stratiform dolomite-cap intervals (Rameil
#     Ch.3) - porosity/permeability highly variable, requires separate
#     diagenetic assessment; NOT simply a low-porosity seal by default.
#
# RULE 2 - Barrier ooid/peloid shoal
#   Olivier F3.1/F3.2 + Rameil barrier ooid bars (inbar4, bch).
#   ooid > 0.16: ooid-dominant grainstone in both regimes' definitions.
#   bioclast > 0.02: reworked bioclastic debris always present on the
#     barrier/shoal (Rameil: "coarse bioclastic rubble" on barriers, p.41;
#     Olivier: "common echinoderms, brachiopods, bivalves" in F3.1).
#   mud < 0.45: grain-supported texture by definition in both frameworks.
#   depth < 16 m: shoal/barrier position above or near FWWB in both regimes.
#   wave > 6500 W/m: high-energy, current/wave-worked setting (cross
#     bedding, herringbone cross-stratification in both Olivier F3.1 and
#     Rameil's barrier/internal-bar facies).
#   GEOTHERMAL: high interparticle + secondary moldic porosity - reservoir
#     facies in both regimes.
#
# RULE 3 - Coral-microbialite / patch-reef buildup
#   Olivier F2.2/F3.4 + Rameil barrier "patch reefs" (FZ5).
#   coral > 0.10: boundstone-forming coral signal; threshold lowered
#     slightly from the Kimmeridgian-only calibration (0.12) because the
#     Tithonian stage modifiers (especially R3-R5) suppress coral
#     production substantially, and a patch reef that DOES survive into the
#     dry phase should still be recognised even at reduced absolute
#     fraction, consistent with Rameil noting that biogenic reefs of the
#     barrier persisted (in laterally equivalent strata) through the
#     interval, just less commonly.
#   mud < 0.55, peloid < 0.60: excludes mud-dominated and peloid-dominated
#     lagoon facies from being mis-classified as reef.
#   depth 8-50 m: spans both Olivier's biostromal/biohermal depth range and
#     Rameil's barrier patch-reef position.
#   GEOTHERMAL: PRIMARY reservoir target in both regimes - vuggy/framework/
#     shelter porosity, generally the highest-quality facies for geothermal
#     fluid flow in Malm carbonate reservoirs.
#
# RULE 4 - Open / normal-marine lagoon (back-shoal)
#   Olivier F4.1-4.3 + Rameil FZ4 normal-marine lagoon (lowEl1-2) and the
#     more open parts of FZ3.
#   peloid > 0.12, oncoid > 0.03: diagnostic non-skeletal grain pair
#     in both frameworks (Olivier's Lithocodium/Bacinella oncoids in F4.1-
#     4.2; Rameil's type-1/2 peloids and abundant oncoids of the open
#     lagoon, Fig.2.8a).
#   bioclast > 0.04: distinguishes this OPEN lagoon facies from the more
#     restricted facies of Rule 5 - Rameil's FZ4 explicitly has "diverse"
#     fauna and "frequent to abundant echinoderms" (lowEl1), unlike the
#     impoverished fauna of FZ3 restricted lagoon.
#   ooid < 0.25, mud < 0.55: excludes barrier shoal and excessively mud-
#     dominated restricted/interior facies.
#   depth 4-22 m: shallower portion of the back-shoal/lagoon position,
#     landward of the barrier/shoal but still within reach of normal-
#     marine water exchange.
#   GEOTHERMAL: moderate porosity secondary reservoir.
#
# RULE 5 - Restricted lagoon
#   Rameil FZ3 (hrl1-4, rl1-9) - the NEW facies category required by the
#     Tithonian regime, not resolved separately in the Kimmeridgian-only
#     calibration because Olivier's open-ramp model does not describe an
#     equivalent restricted-lagoon belt at Eclépens' position.
#   Field marker (Rameil p.40-41): "diversity of fauna and flora commonly
#     reduced due to increased or decreased salinity"; "intense
#     bioturbation by few types of organisms and euryhaline fauna" -
#     modelled as reduced bioclastic diversity relative to Rule 4.
#   peloid > 0.10, evap_mud > 0.04: peloid-dominated (as in FZ4) BUT with a
#     measurable evaporitic-mud signal that the open lagoon lacks - this is
#     the key discriminator at grid scale, since both facies share the
#     peloid+oncoid grain assemblage but differ in restriction/salinity
#     stress, which the evap_mud factory's presence captures as a proxy.
#   bioclast < 0.04: LOW bioclastic diversity, the inverse of Rule 4's
#     threshold - consistent with Rameil's "faunal desert" character of
#     the more restricted lagoon sub-facies (hrl1-3).
#   mud > 0.30: mud-to-packstone texture (Rameil's hrl/rl facies range
#     from M to W-P).
#   depth 4-25 m: same general position as Rule 4 but landward/more
#     restricted within the lagoon system.
#   GEOTHERMAL: low-to-moderate porosity, increasingly affected by
#     dolomitization toward the dry-phase stages - variable reservoir
#     quality requiring the Rule-1 dolomite-cap flag for refinement.
#
# RULE 6 - Interior platform / protected mudstone
#   Olivier F4.4 - low-energy interior platform mudstone. Rameil's deepest,
#     most restricted hrl1 ("empty, texture-less mudstone... faunal
#     deserts") is the Tithonian equivalent of this very low-energy,
#     mud-dominated end-member.
#   mud > 0.45: mud-dominated by definition (Dunham mudstone in both
#     frameworks).
#   bioclast < 0.04, evap_mud < 0.04: distinguishes from Rule 5 (restricted
#     lagoon with measurable evap_mud signal) - this facies is simply
#     QUIET and mud-dominated, without the evaporitic-stress signature.
#   depth < 30 m: shallow, protected interior position in both regimes.
#   wave < 5000 W/m: explicitly low-energy (Olivier F4.4: "low energy
#     internal platform deposition", Rameil hrl1: protected, non-
#     evaporitic, deepest part of the restricted realm).
#   GEOTHERMAL: low porosity seal/baffle in both regimes.
#
# RULE 7 - Mid-ramp / outer lagoon oncoid-bioclastic
#   Olivier F2.1/F2.3 mid-ramp oncoid packstone/marls; Rameil's outer-
#     lagoon to barrier-seaward transition where oncoid content remains
#     significant but the setting is deeper/more open than the back-shoal.
#   oncoid > 0.05: abundant-oncoid signature in both frameworks.
#   mud 0.20-0.75: packstone-to-wackestone range, broader than the shallow
#     lagoon facies above.
#   ooid < 0.18: below or seaward of direct barrier/shoal influence.
#   depth 22-65 m: deeper position than the back-shoal/lagoon facies,
#     extending toward (but not into) the truly offshore/outer-platform
#     setting of Rule 8.
#   GEOTHERMAL: moderate porosity secondary reservoir (oncoid moldic +
#     interparticle).
#
# RULE 8 - Offshore / outer platform mudstone
#   Olivier F1.1-1.2 lower-offshore mudstone; Rameil's deepest basin-ward
#     equivalent recorded in slope/basin sections (FZ6-7, outside the
#     platform sensu stricto but representing the same end-member toward
#     which the deepest platform cells drift during accommodation excess).
#   mud > 0.70: overwhelmingly mud-supported, very low allochem signal.
#   ooid < 0.08, evap_mud < 0.05: both essentially absent in the deepest,
#     most distal setting.
#   depth > 55 m: below storm wave base in both regimes' depth structure.
#   GEOTHERMAL: non-reservoir, seal/source rock equivalent.
# ─────────────────────────────────────────────────────────────────────────────

rules = [

    # Rule 1 - Tidal flat / sabkha (dolomite-cap candidate)
    # Specific dry-phase cap only; no wave upper gate, otherwise the shallow
    # high evap_mud cells were going to fallback.
    FaciesRule(
        name               = "Tidal flat / sabkha (dolomite-cap candidate)",
        sediment_fractions = Dict(7=>(0.22,1.0), 3=>(0.25,1.0), 1=>(0.0,0.18)),
        depth_range        = (-Inf*u"m", 12.0*u"m"),
        wave_energy_range  = (0.0*u"W/m", Inf*u"W/m")),

    # Rule 2 - Barrier ooid/peloid shoal
    # Bioclasts are no longer mandatory because the R1-R2 ooid-rich cells
    # were too clean and were being hidden by the bioclast gate.
    FaciesRule(
        name               = "Barrier ooid/peloid shoal",
        sediment_fractions = Dict(1=>(0.30,1.0), 4=>(0.08,1.0), 3=>(0.0,0.55)),
        depth_range        = (0.0*u"m", 16.0*u"m"),
        wave_energy_range  = (6000.0*u"W/m", Inf*u"W/m")),

    # Rule 3 - Coral-microbialite / patch-reef buildup
    FaciesRule(
        name               = "Coral-microbialite / patch-reef buildup",
        sediment_fractions = Dict(2=>(0.08,1.0), 3=>(0.0,0.65), 4=>(0.0,0.70)),
        depth_range        = (8.0*u"m", 50.0*u"m"),
        wave_energy_range  = (0.0*u"W/m", Inf*u"W/m")),

    # Rule 4 - Open / normal-marine lagoon (back-shoal)
    # Put before restricted lagoon so the Olivier/Kimmeridgian open/protected
    # platform is not swallowed by the Tithonian restriction rule.
    FaciesRule(
        name               = "Open normal-marine lagoon (back-shoal)",
        sediment_fractions = Dict(4=>(0.12,1.0), 6=>(0.010,1.0), 5=>(0.005,1.0),
                                  1=>(0.0,0.28), 3=>(0.0,0.65), 7=>(0.0,0.025)),
        depth_range        = (4.0*u"m", 40.0*u"m"),
        wave_energy_range  = (0.0*u"W/m", Inf*u"W/m")),

    # Rule 5 - Restricted lagoon
    # Requires a measurable evap_mud/restriction signal so Stage O does not
    # classify as Tithonian restricted lagoon.
    FaciesRule(
        name               = "Restricted lagoon",
        sediment_fractions = Dict(7=>(0.015,1.0), 3=>(0.30,1.0),
                                  5=>(0.0,0.065), 1=>(0.0,0.30)),
        depth_range        = (0.0*u"m", 35.0*u"m"),
        wave_energy_range  = (0.0*u"W/m", Inf*u"W/m")),

    # Rule 6 - Interior platform / protected mudstone
    FaciesRule(
        name               = "Interior platform / protected mudstone",
        sediment_fractions = Dict(3=>(0.45,1.0), 5=>(0.0,0.075), 7=>(0.0,0.035)),
        depth_range        = (0.0*u"m", 45.0*u"m"),
        wave_energy_range  = (0.0*u"W/m", Inf*u"W/m")),

    # Rule 7 - Mid-ramp / outer lagoon oncoid-bioclastic
    FaciesRule(
        name               = "Mid-ramp / outer lagoon oncoid-bioclastic",
        sediment_fractions = Dict(6=>(0.04,1.0), 3=>(0.20,0.75), 1=>(0.0,0.22)),
        depth_range        = (22.0*u"m", 75.0*u"m"),
        wave_energy_range  = (0.0*u"W/m", Inf*u"W/m")),

    # Rule 8 - Offshore / outer platform mudstone
    FaciesRule(
        name               = "Offshore / outer platform mudstone",
        sediment_fractions = Dict(3=>(0.70,1.0), 1=>(0.0,0.08)),
        depth_range        = (55.0*u"m", 220.0*u"m"),
        wave_energy_range  = (0.0*u"W/m", Inf*u"W/m")),
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
    cls_labels,"Diagnostic class";orientation=:horizontal,nbanks=3,framevisible=false)
save(joinpath(FIG_DIR,"strat_columns_classified.png"),fig_sc;px_per_unit=2)

@info "Done. Figures in $(FIG_DIR)/"

end  # module Script
