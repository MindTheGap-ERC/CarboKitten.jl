# export_carbokitten_to_petrel_grdecl.jl
# Fixed version: 7 production factories + 9 classified facies for final Eclépens calibration
#
# Purpose:
#   Read CarboKitten HDF5 outputs and export one Petrel/ECLIPSE-style GRDECL
#   corner-point grid containing:
#     - geometry from a GOCAD TSurf top surface
#     - one layer per saved CarboKitten output interval
#     - dominant production factory property
#     - dominant classified facies property
#     - optional production/classification fraction properties
#
# Expected HDF5 structure:
#   production HDF5:
#     topography/production          (nt, ny, nx, nprod)
#     topography/sediment_thickness  (nt, ny, nx)
#
#   classified HDF5:
#     classified/production          (nt, ny, nx, nclass)
#     classified/sediment_thickness  (nt, ny, nx)
#
# Usage from repository root:
#   julia --project=. examples/model/export_carbokitten_to_petrel_grdecl.jl

using HDF5
using Printf
using Statistics
using LinearAlgebra

# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------

const PROD_H5     = "data/output/eclepens_withoutca.h5"
const CLASS_H5    = "data/output/eclepens_withoutca_classified.h5"
const TSURF_FILE  = "data/input/coordinates_malm.ts"
const OUT_GRDECL  = "data/output/eclepens_carbokitten_petrel_grid_keywords_7prod_9facies_poro_claytight_above_lmalm_shiftup265.GRDECL"

# Use :bbox for a non-rotated Petrel grid covering the TSurf extent.
# Use :pca if your Petrel grid is rotated along the main TSurf trend.
const GRID_MODE = :bbox       # :bbox or :pca

# The uploaded TSurf / bathymetry surface is the Lower Malm reference surface.
# The CarboKitten package is deposited ABOVE this surface, not below it.
# The uploaded TSurf says ZPOSITIVE Elevation and contains negative Z values.
# Use :elevation to preserve those Petrel elevations.
# If Petrel imports the grid upside down, switch to :depth.
const Z_MODE = :depth         # positive depth for GRDECL/Petrel

# Shift the whole exported GRDECL grid upward by this amount.
# In depth mode, upward means subtracting from depth.
# In elevation mode, upward means adding to elevation.
const VERTICAL_SHIFT_UP_M = 265.0

# HDF5 output frame 1 is normally zero thickness. Keep false for Petrel export.
const INCLUDE_ZERO_FIRST_FRAME = false

# Write continuous fraction properties as well as dominant integer properties.
# This makes the GRDECL larger but much more useful in Petrel.
const WRITE_FRACTIONS = true

# Cells thinner than this are inactive.
const MIN_CELL_THICKNESS_M = 0.05

# Extra vertical length added to COORD pillars below the deepest ZCORN.
const PILLAR_EXTRA_M = 50.0

# Invert the HDF5 X/Y indices when mapping to Petrel coordinates.
# For Petrel/ECLIPSE left-handed I/J orientation, reverse one horizontal axis only.
const FLIP_X = false
const FLIP_Y = true

# IDW interpolation settings for the TSurf top surface.
const IDW_K = 12
const IDW_POWER = 2.0

# Property labels.
# These must match the final calibrated run-ecl.jl rule order.
const PROD_NAMES = [
    "ooids",
    "corals",
    "mud",
    "peloids",
    "bioclasts_intraclasts",
    "oncoids",
    "evap_mud",
]

const CLASS_NAMES = [
    "Tidal flat / sabkha (dolomite-cap candidate)",
    "Barrier ooid/peloid shoal",
    "Coral-microbialite / patch-reef buildup",
    "Open normal-marine lagoon (back-shoal)",
    "Restricted lagoon",
    "Interior platform / protected mudstone",
    "Mid-ramp / outer lagoon oncoid-bioclastic",
    "Offshore / outer platform mudstone",
    "fallback",
]

const PROD_DOM_KEY      = "LITHO"
const CLASS_DOM_KEY     = "FACIES"
const PROD_ALIAS_KEY    = "FIPNUM"
const CLASS_ALIAS_KEY   = "SATNUM"
const WRITE_CLASSIC_ALIASES = true

# Porosity-quality class property.
# PORO is exported as integer classes because Petrel imports this keyword reliably.
# Treat it as a categorical property/template after import, not as calibrated porosity:
#   1 = CLAY     clay-rich/marly or very mud-rich tight carbonate proxy
#   2 = TIGHTLS  tight limestone
#   3 = MEDPOR   medium porous limestone
#   4 = HIGHPOR  highly porous limestone
const PORO_KEY = "PORO"
const PORO_NAMES = ["CLAY", "TIGHTLS", "MEDPOR", "HIGHPOR"]

# Keep GRDECL property keywords short and Petrel/ECLIPSE-friendly.
const PROD_FRAC_KEYS = ["OOID", "CORAL", "MUD", "PELOID", "BIOCL", "ONCOID", "EVPMUD"]
const CLS_FRAC_KEYS  = ["SABKHA", "BSHOAL", "REEF", "OPENLAG", "RESTLAG",
                        "INMUD", "MIDRAMP", "OFFSH", "UNCLASS"]

# ---------------------------------------------------------------------------
# TSurf parsing and surface interpolation
# ---------------------------------------------------------------------------

function read_tsurf_vertices(path::AbstractString)
    xs = Float64[]
    ys = Float64[]
    zs = Float64[]

    open(path, "r") do io
        for line in eachline(io)
            s = split(strip(line))
            if length(s) >= 5 && (s[1] == "VRTX" || s[1] == "PVRTX")
                push!(xs, parse(Float64, s[3]))
                push!(ys, parse(Float64, s[4]))
                push!(zs, parse(Float64, s[5]))
            end
        end
    end

    if isempty(xs)
        error("No VRTX/PVRTX vertices found in $path")
    end

    return hcat(xs, ys, zs)
end

function make_pillar_xy(vertices::Matrix{Float64}, nx::Int, ny::Int; mode::Symbol=:bbox)
    xy = vertices[:, 1:2]

    if mode == :bbox
        xmin, xmax = extrema(xy[:, 1])
        ymin, ymax = extrema(xy[:, 2])

        xs = collect(range(xmin, xmax; length=nx+1))
        ys = collect(range(ymin, ymax; length=ny+1))

        X = zeros(Float64, ny+1, nx+1)
        Y = zeros(Float64, ny+1, nx+1)

        for j in 1:ny+1, i in 1:nx+1
            ii = FLIP_X ? (nx + 2 - i) : i
            jj = FLIP_Y ? (ny + 2 - j) : j
            X[j, i] = xs[ii]
            Y[j, i] = ys[jj]
        end

        return X, Y

    elseif mode == :pca
        center = vec(mean(xy; dims=1))
        A = xy .- center'
        C = cov(A; dims=1)
        eig = eigen(Symmetric(C))
        order = sortperm(eig.values; rev=true)
        V = eig.vectors[:, order]

        uv = A * V
        umin, umax = extrema(uv[:, 1])
        vmin, vmax = extrema(uv[:, 2])

        us = collect(range(umin, umax; length=nx+1))
        vs = collect(range(vmin, vmax; length=ny+1))

        X = zeros(Float64, ny+1, nx+1)
        Y = zeros(Float64, ny+1, nx+1)

        for j in 1:ny+1, i in 1:nx+1
            ii = FLIP_X ? (nx + 2 - i) : i
            jj = FLIP_Y ? (ny + 2 - j) : j
            p = center .+ V * [us[ii], vs[jj]]
            X[j, i] = p[1]
            Y[j, i] = p[2]
        end

        return X, Y

    else
        error("Unknown GRID_MODE=$mode. Use :bbox or :pca.")
    end
end

function build_spatial_index(vertices::Matrix{Float64}, cellsize::Float64)
    xmin = minimum(vertices[:, 1])
    ymin = minimum(vertices[:, 2])

    bins = Dict{Tuple{Int, Int}, Vector{Int}}()

    for n in axes(vertices, 1)
        ix = floor(Int, (vertices[n, 1] - xmin) / cellsize)
        iy = floor(Int, (vertices[n, 2] - ymin) / cellsize)
        key = (ix, iy)
        if !haskey(bins, key)
            bins[key] = Int[]
        end
        push!(bins[key], n)
    end

    return bins, xmin, ymin
end

function idw_z(x::Float64, y::Float64, vertices::Matrix{Float64},
               bins::Dict{Tuple{Int, Int}, Vector{Int}},
               xmin::Float64, ymin::Float64, cellsize::Float64;
               k::Int=IDW_K, power::Float64=IDW_POWER)

    ix = floor(Int, (x - xmin) / cellsize)
    iy = floor(Int, (y - ymin) / cellsize)

    cand = Int[]
    maxring = 50

    for ring in 0:maxring
        empty!(cand)
        for bx in (ix-ring):(ix+ring), by in (iy-ring):(iy+ring)
            append!(cand, get(bins, (bx, by), Int[]))
        end
        if length(cand) >= k || ring == maxring
            break
        end
    end

    if isempty(cand)
        error("No TSurf interpolation candidates near x=$x, y=$y")
    end

    d2 = [(vertices[n, 1] - x)^2 + (vertices[n, 2] - y)^2 for n in cand]

    nearest_local = partialsortperm(d2, 1:min(k, length(d2)))

    num = 0.0
    den = 0.0

    for iloc in nearest_local
        n = cand[iloc]
        dist = sqrt(d2[iloc])
        if dist < 1e-9
            return vertices[n, 3]
        end
        w = 1.0 / dist^power
        num += w * vertices[n, 3]
        den += w
    end

    return num / den
end

function interpolate_top_surface(vertices::Matrix{Float64}, X::Matrix{Float64}, Y::Matrix{Float64})
    ny1, nx1 = size(X)

    # Bin size roughly adapted to Petrel grid resolution.
    dx = (maximum(X) - minimum(X)) / max(nx1 - 1, 1)
    dy = (maximum(Y) - minimum(Y)) / max(ny1 - 1, 1)
    cellsize = max(dx, dy) * 4.0

    bins, xmin, ymin = build_spatial_index(vertices, cellsize)

    Z = zeros(Float64, ny1, nx1)
    for j in 1:ny1, i in 1:nx1
        Z[j, i] = idw_z(X[j, i], Y[j, i], vertices, bins, xmin, ymin, cellsize)
    end

    return Z
end

# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

function centers_to_corners(A::Matrix{Float64})
    ny, nx = size(A)
    C = zeros(Float64, ny+1, nx+1)
    W = zeros(Float64, ny+1, nx+1)

    for j in 1:ny, i in 1:nx
        v = A[j, i]
        C[j,   i]   += v; W[j,   i]   += 1.0
        C[j,   i+1] += v; W[j,   i+1] += 1.0
        C[j+1, i]   += v; W[j+1, i]   += 1.0
        C[j+1, i+1] += v; W[j+1, i+1] += 1.0
    end

    return C ./ W
end

function z_at(lower_malm_z::Float64, cumulative_thickness::Float64)
    # lower_malm_z is the Petrel/TSurf elevation of the Lower Malm reference surface.
    # Positive cumulative_thickness is deposited ABOVE Lower Malm.
    #
    # The final VERTICAL_SHIFT_UP_M moves the whole exported grid upward.
    #
    # If Z_MODE = :depth and lower_malm_z = -800 m:
    #   cumulative_thickness = 0 m, shift = 265 m  -> z = 535 m depth
    #   cumulative_thickness = 20 m, shift = 265 m -> z = 515 m depth
    if Z_MODE == :elevation
        return lower_malm_z + cumulative_thickness + VERTICAL_SHIFT_UP_M
    elseif Z_MODE == :depth
        return -lower_malm_z - cumulative_thickness - VERTICAL_SHIFT_UP_M
    else
        error("Unknown Z_MODE=$Z_MODE. Use :elevation or :depth.")
    end
end

function write_values(io, values; perline::Int=6, fmt::String="%.6f")
    # @printf requires a literal format string in Julia 1.12.
    # Use Printf.Format for runtime-selected formats.
    fmtobj = Printf.Format(" " * fmt)

    n = 0
    for v in values
        print(io, Printf.format(fmtobj, v))
        n += 1
        if n % perline == 0
            print(io, "\n")
        end
    end
    if n % perline != 0
        print(io, "\n")
    end
    print(io, "/\n")
end

function write_int_values(io, values; perline::Int=20)
    n = 0
    for v in values
        @printf(io, " %d", v)
        n += 1
        if n % perline == 0
            print(io, "\n")
        end
    end
    if n % perline != 0
        print(io, "\n")
    end
    print(io, "/\n")
end

function dominant_index(v::AbstractVector{<:Real})
    s = sum(v)
    if !(s > 0.0)
        return 0
    end
    imax = 1
    vmax = v[1]
    for i in 2:length(v)
        if v[i] > vmax
            vmax = v[i]
            imax = i
        end
    end
    return imax
end

function fraction_value(v::AbstractVector{<:Real}, idx::Int)
    s = sum(v)
    return s > 0.0 ? Float64(v[idx] / s) : 0.0
end

function poro_class(cls_vec::AbstractVector{<:Real}, prod_vec::AbstractVector{<:Real})
    cls = dominant_index(cls_vec)
    cls == 0 && return 0

    s = sum(prod_vec)
    if !(s > 0.0)
        # Conservative fallback from depositional facies only.
        return cls == 2 ? 3 :   # barrier shoal -> MEDPOR unless clean-grain override says HIGHPOR
               cls == 3 ? 3 :   # coral patch reef -> MEDPOR unless clean-grain override says HIGHPOR
               cls == 4 ? 2 :   # open/back-shoal lagoon -> TIGHTLS
               cls == 7 ? 2 :   # mid-ramp/oncoid-bioclastic -> TIGHTLS
               cls == 8 ? 1 :   # offshore mudstone -> CLAY
               cls == 6 ? 1 :   # protected mudstone -> CLAY
               2                # all others -> TIGHTLS
    end

    ooid    = length(prod_vec) >= 1 ? prod_vec[1] / s : 0.0
    coral   = length(prod_vec) >= 2 ? prod_vec[2] / s : 0.0
    mud     = length(prod_vec) >= 3 ? prod_vec[3] / s : 0.0
    peloid  = length(prod_vec) >= 4 ? prod_vec[4] / s : 0.0
    bioc    = length(prod_vec) >= 5 ? prod_vec[5] / s : 0.0
    oncoid  = length(prod_vec) >= 6 ? prod_vec[6] / s : 0.0
    evap    = length(prod_vec) >= 7 ? prod_vec[7] / s : 0.0
    grain   = ooid + coral + peloid + bioc + oncoid
    clean_grain = ooid + coral + bioc

    # Conservative Eclepens Malm well-oriented mapping:
    # mostly CLAY and TIGHTLS; MEDPOR/HIGHPOR only where composition is clearly favourable.
    #
    # CLAY = clay-rich/marly or very mud-rich tight carbonate proxy.
    if cls == 8
        return 1
    elseif cls == 6 && mud >= 0.50
        return 1
    elseif mud >= 0.62 && grain <= 0.38
        return 1
    elseif cls == 5 && mud >= 0.55 && evap < 0.15
        return 1
    end

    # HIGHPOR only for clean shoal/reef-prone cells.
    if (cls == 2 || cls == 3) && clean_grain >= 0.32 && mud <= 0.45
        return 4
    end

    # MEDPOR only for genuinely grainier or dolomitization-prone intervals.
    if evap >= 0.15 && mud <= 0.70
        return 3
    elseif (cls == 2 || cls == 3) && grain >= 0.35 && mud <= 0.60
        return 3
    elseif (cls == 4 || cls == 7) && grain >= 0.42 && mud <= 0.52
        return 3
    end

    # Everything else is tight limestone by default.
    return 2
end

function read_grid_lengths(h5path::AbstractString)
    h5open(h5path, "r") do f
        if haskey(f, "input/x") && haskey(f, "input/y")
            nx = length(read(f["input/x"]))
            ny = length(read(f["input/y"]))
            return nx, ny
        end
    end
    return nothing
end

function normalize_cumulative(raw::AbstractArray, nx::Int, ny::Int)
    sz = size(raw)

    if length(sz) != 3
        error("sediment_thickness must be 3D, got size $(sz)")
    end

    # Target internal order: (nt, ny, nx)
    if sz[2] == ny && sz[3] == nx
        # Already (nt, ny, nx)
        return Float64.(raw)
    elseif sz[1] == nx && sz[2] == ny
        # Julia/HDF5 order from CarboKitten: (nx, ny, nt)
        return Float64.(permutedims(raw, (3, 2, 1)))
    elseif sz[1] == ny && sz[2] == nx
        # Alternate order: (ny, nx, nt)
        return Float64.(permutedims(raw, (3, 1, 2)))
    elseif sz[1] == nx && sz[3] == ny
        # Alternate order: (nx, nt, ny)
        return Float64.(permutedims(raw, (2, 3, 1)))
    elseif sz[1] == ny && sz[3] == nx
        # Alternate order: (ny, nt, nx)
        return Float64.(permutedims(raw, (2, 1, 3)))
    else
        error("Could not normalize sediment_thickness size $(sz) with nx=$nx ny=$ny")
    end
end

function normalize_component_data(raw::AbstractArray, nx::Int, ny::Int)
    sz = size(raw)

    if length(sz) != 4
        error("production/classified data must be 4D, got size $(sz)")
    end

    # Target internal order: (nt, ny, nx, nf)
    if sz[2] == ny && sz[3] == nx
        # Already (nt, ny, nx, nf)
        return Float64.(raw)
    elseif sz[2] == nx && sz[3] == ny
        # Julia/HDF5 order from CarboKitten: (nf, nx, ny, nt)
        return Float64.(permutedims(raw, (4, 3, 2, 1)))
    elseif sz[1] == nx && sz[2] == ny
        # Alternate order: (nx, ny, nf, nt)
        return Float64.(permutedims(raw, (4, 2, 1, 3)))
    elseif sz[1] == ny && sz[2] == nx
        # Alternate order: (ny, nx, nf, nt)
        return Float64.(permutedims(raw, (4, 1, 2, 3)))
    elseif sz[1] == nx && sz[3] == ny
        # Alternate order: (nx, nf, ny, nt)
        return Float64.(permutedims(raw, (4, 3, 1, 2)))
    elseif sz[1] == ny && sz[3] == nx
        # Alternate order: (ny, nf, nx, nt)
        return Float64.(permutedims(raw, (4, 1, 3, 2)))
    else
        error("Could not normalize component data size $(sz) with nx=$nx ny=$ny")
    end
end

function infer_grid_lengths(prod_raw::AbstractArray, cum_raw::AbstractArray)
    # Prefer input/x and input/y, but this fallback handles older files.
    csz = size(cum_raw)
    psz = size(prod_raw)

    if length(csz) == 3
        # Common Julia/CarboKitten order: cumulative = (nx, ny, nt)
        # Common Python-visible order: cumulative = (nt, ny, nx)
        # Pick the two largest non-time-like dimensions.
        dims = collect(csz)
        sorted = sort(dims; rev=true)
        ny = sorted[1]
        nx = sorted[2]
        return nx, ny
    end

    error("Cannot infer grid size from cumulative size $(csz) and production size $(psz)")
end

# ---------------------------------------------------------------------------
# GRDECL writer
# ---------------------------------------------------------------------------

function build_zcorn(top_surface::Matrix{Float64}, cumulative::Array{Float64, 3},
                     kframes::Vector{Int})
    # cumulative has dimensions nt, ny, nx.
    #
    # The Lower Malm surface is the BASE of the modelled carbonate package.
    # Therefore each depositional interval lies ABOVE the Lower Malm surface:
    #
    #   base of interval = cumulative thickness at previous frame
    #   top  of interval = cumulative thickness at current frame
    #
    # kframes is ordered top-to-base before export, so Petrel K=1 is the youngest/top layer.
    nt, ny, nx = size(cumulative)
    nz = length(kframes)

    zcorn = zeros(Float64, 2*nx, 2*ny, 2*nz)
    zero_corners = centers_to_corners(zeros(Float64, ny, nx))

    for kk in 1:nz
        kt = kframes[kk]

        top_corners = centers_to_corners(cumulative[kt, :, :])
        base_corners = kt > 1 ? centers_to_corners(cumulative[kt - 1, :, :]) : zero_corners

        for j in 1:ny, i in 1:nx
            top_sw = z_at(top_surface[j,   i],   top_corners[j,   i])
            top_se = z_at(top_surface[j,   i+1], top_corners[j,   i+1])
            top_nw = z_at(top_surface[j+1, i],   top_corners[j+1, i])
            top_ne = z_at(top_surface[j+1, i+1], top_corners[j+1, i+1])

            bot_sw = z_at(top_surface[j,   i],   base_corners[j,   i])
            bot_se = z_at(top_surface[j,   i+1], base_corners[j,   i+1])
            bot_nw = z_at(top_surface[j+1, i],   base_corners[j+1, i])
            bot_ne = z_at(top_surface[j+1, i+1], base_corners[j+1, i+1])

            zcorn[2*i-1, 2*j-1, 2*kk-1] = top_sw
            zcorn[2*i,   2*j-1, 2*kk-1] = top_se
            zcorn[2*i-1, 2*j,   2*kk-1] = top_nw
            zcorn[2*i,   2*j,   2*kk-1] = top_ne

            zcorn[2*i-1, 2*j-1, 2*kk] = bot_sw
            zcorn[2*i,   2*j-1, 2*kk] = bot_se
            zcorn[2*i-1, 2*j,   2*kk] = bot_nw
            zcorn[2*i,   2*j,   2*kk] = bot_ne
        end
    end

    return vec(zcorn)
end

function build_coord(X::Matrix{Float64}, Y::Matrix{Float64}, top_surface::Matrix{Float64},
                     cumulative::Array{Float64, 3})
    ny1, nx1 = size(X)
    nt, ny, nx = size(cumulative)

    # Top of the exported grid is the final cumulative depositional top.
    # Base of the exported grid is the Lower Malm surface, slightly extended down the pillars.
    max_cum_center = cumulative[end, :, :]
    max_cum_corner = centers_to_corners(max_cum_center)

    coord = Float64[]

    for j in 1:ny1, i in 1:nx1
        ztop = z_at(top_surface[j, i], max_cum_corner[j, i])
        zbot = z_at(top_surface[j, i], -PILLAR_EXTRA_M)
        append!(coord, (X[j, i], Y[j, i], ztop, X[j, i], Y[j, i], zbot))
    end

    return coord
end

function build_actnum(cumulative::Array{Float64, 3}, kframes::Vector{Int})
    nt, ny, nx = size(cumulative)
    vals = Int[]
    zero_surface = zeros(Float64, ny, nx)

    for kt in kframes
        prev = kt > 1 ? cumulative[kt - 1, :, :] : zero_surface
        cur = cumulative[kt, :, :]
        thick = cur .- prev

        for j in 1:ny, i in 1:nx
            push!(vals, thick[j, i] > MIN_CELL_THICKNESS_M ? 1 : 0)
        end
    end

    return vals
end

function build_dominant_property(data::Array{Float64, 4}, kframes::Vector{Int})
    nt, ny, nx, nf = size(data)
    vals = Int[]

    for kt in kframes
        for j in 1:ny, i in 1:nx
            push!(vals, dominant_index(@view data[kt, j, i, :]))
        end
    end

    return vals
end

function build_fraction_property(data::Array{Float64, 4}, kframes::Vector{Int}, ifac::Int)
    nt, ny, nx, nf = size(data)
    vals = Float64[]

    for kt in kframes
        for j in 1:ny, i in 1:nx
            push!(vals, fraction_value(@view(data[kt, j, i, :]), ifac))
        end
    end

    return vals
end

function build_poro_property(cls::Array{Float64, 4}, prod::Array{Float64, 4}, kframes::Vector{Int})
    nt, ny, nx, ncls = size(cls)
    ntp, nyp, nxp, nprod = size(prod)

    if (ntp, nyp, nxp) != (nt, ny, nx)
        error("Cannot build PETRO: classified grid $(size(cls)) and production grid $(size(prod)) do not match.")
    end

    vals = Int[]

    for kt in kframes
        for j in 1:ny, i in 1:nx
            push!(vals, poro_class(@view(cls[kt, j, i, :]), @view(prod[kt, j, i, :])))
        end
    end

    return vals
end

function export_grdecl()
    println("Reading HDF5...")
    prod_raw = h5open(PROD_H5, "r") do f
        read(f["topography/production"])
    end

    cum_raw = h5open(PROD_H5, "r") do f
        read(f["topography/sediment_thickness"])
    end

    cls_raw = h5open(CLASS_H5, "r") do f
        read(f["classified/production"])
    end

    grid_lengths = read_grid_lengths(PROD_H5)
    nx, ny = grid_lengths === nothing ? infer_grid_lengths(prod_raw, cum_raw) : grid_lengths

    println("Raw dataset sizes:")
    println("  topography/production         = $(size(prod_raw))")
    println("  topography/sediment_thickness = $(size(cum_raw))")
    println("  classified/production         = $(size(cls_raw))")
    println("  input grid                    = nx=$nx ny=$ny")

    prod = normalize_component_data(prod_raw, nx, ny)
    cum  = normalize_cumulative(cum_raw, nx, ny)
    cls  = normalize_component_data(cls_raw, nx, ny)

    nt, ny_norm, nx_norm, nprod = size(prod)
    ntc, nyc, nxc, ncls = size(cls)

    if size(cum) != (nt, ny_norm, nx_norm)
        error("Normalized sediment_thickness size $(size(cum)) does not match normalized production size $(size(prod)).")
    end
    if (ntc, nyc, nxc) != (nt, ny_norm, nx_norm)
        error("Normalized classified/production size $(size(cls)) does not match normalized production grid $(size(prod)).")
    end

    nx = nx_norm
    ny = ny_norm

    interval_frames = INCLUDE_ZERO_FIRST_FRAME ? collect(1:nt) : collect(2:nt)

    # Petrel expects K layers ordered from top to base.
    # Because deposition is above the Lower Malm base, the final frame is the topmost layer.
    kframes = reverse(interval_frames)
    nz = length(kframes)

    println("Normalized internal sizes:")
    println("  production = $(size(prod))  interpreted as (nt, ny, nx, nprod)")
    println("  cumulative = $(size(cum))   interpreted as (nt, ny, nx)")
    println("  classified = $(size(cls))   interpreted as (nt, ny, nx, nclass)")
    println("Grid dimensions exported: nx=$nx ny=$ny nz=$nz")
    println("Production factories: $nprod")
    println("Classified classes: $ncls")

    if nprod != length(PROD_NAMES)
        error("HDF5 production has $nprod factories but PROD_NAMES has $(length(PROD_NAMES)). Update PROD_NAMES/PROD_FRAC_KEYS.")
    end
    if ncls != length(CLASS_NAMES)
        error("HDF5 classified data has $ncls classes but CLASS_NAMES has $(length(CLASS_NAMES)). Update CLASS_NAMES/CLS_FRAC_KEYS.")
    end

    println("Reading and interpolating TSurf geometry...")
    verts = read_tsurf_vertices(TSURF_FILE)
    X, Y = make_pillar_xy(verts, nx, ny; mode=GRID_MODE)
    top_surface = interpolate_top_surface(verts, X, Y)

    println("Building GRDECL arrays...")
    coord = build_coord(X, Y, top_surface, cum)
    zcorn = build_zcorn(top_surface, cum, kframes)
    actnum = build_actnum(cum, kframes)
    prod_dom = build_dominant_property(prod, kframes)
    cls_dom = build_dominant_property(cls, kframes)
    poro = build_poro_property(cls, prod, kframes)

    println("Writing $OUT_GRDECL ...")
    mkpath(dirname(OUT_GRDECL))

    open(OUT_GRDECL, "w") do io
        println(io, "-- Exported from CarboKitten HDF5")
        println(io, "-- Production HDF5: $PROD_H5")
        println(io, "-- Classified HDF5: $CLASS_H5")
        println(io, "-- TSurf geometry: $TSURF_FILE")
        println(io, "-- GRID_MODE=$GRID_MODE Z_MODE=$Z_MODE")
        println(io, "-- Lower Malm is treated as the base of the modelled package.")
        println(io, "-- Vertical shift upward applied: $(VERTICAL_SHIFT_UP_M) m")
        println(io, "-- Depositional thickness is exported upward above Lower Malm, not downward below it.")
        println(io, "-- I/J handedness fix: FLIP_X=$FLIP_X FLIP_Y=$FLIP_Y")
        println(io, "-- Single-axis flip changes handedness; double-axis flip only rotates 180 degrees.")
        println(io, "-- Grid keywords: SPECGRID, COORD, ZCORN, ACTNUM")
        println(io, "-- Display properties: FACIES, LITHO, PORO")
        println(io, "-- Classic aliases: SATNUM = FACIES, FIPNUM = LITHO")
        println(io, "--")
        println(io, "-- LITHO/FIPNUM values:")
        for (i, name) in enumerate(PROD_NAMES)
            println(io, "--   $i = $name")
        end
        println(io, "--")
        println(io, "-- FACIES/SATNUM values:")
        for (i, name) in enumerate(CLASS_NAMES)
            println(io, "--   $i = $name")
        end
        println(io, "--   0 = inactive/no deposit")
        println(io, "--")
        println(io, "-- PORO values:")
        for (i, name) in enumerate(PORO_NAMES)
            println(io, "--   $i = $name")
        end
        println(io, "--   0 = inactive/no deposit")
        println(io)

        println(io, "SPECGRID")
        @printf(io, " %d %d %d 1 F /\n\n", nx, ny, nz)

        println(io, "COORD")
        write_values(io, coord; perline=6, fmt="%.6f")
        println(io)

        println(io, "ZCORN")
        write_values(io, zcorn; perline=8, fmt="%.6f")
        println(io)

        println(io, "ACTNUM")
        write_int_values(io, actnum; perline=30)
        println(io)

        # Meaningful Petrel-style discrete properties.
        println(io, CLASS_DOM_KEY)
        write_int_values(io, cls_dom; perline=30)
        println(io)

        println(io, PROD_DOM_KEY)
        write_int_values(io, prod_dom; perline=30)
        println(io)

        println(io, PORO_KEY)
        write_int_values(io, poro; perline=30)
        println(io)

        # Classic ECLIPSE aliases for robust Petrel import/display.
        if WRITE_CLASSIC_ALIASES
            println(io, CLASS_ALIAS_KEY)
            write_int_values(io, cls_dom; perline=30)
            println(io)

            println(io, PROD_ALIAS_KEY)
            write_int_values(io, prod_dom; perline=30)
            println(io)
        end

        if WRITE_FRACTIONS
            for ifac in 1:nprod
                key = ifac <= length(PROD_FRAC_KEYS) ? PROD_FRAC_KEYS[ifac] : "P_F$(ifac)"
                println(io, key)
                write_values(io, build_fraction_property(prod, kframes, ifac); perline=8, fmt="%.5f")
                println(io)
            end

            for icls in 1:ncls
                key = icls <= length(CLS_FRAC_KEYS) ? CLS_FRAC_KEYS[icls] : "C_F$(icls)"
                println(io, key)
                write_values(io, build_fraction_property(cls, kframes, icls); perline=8, fmt="%.5f")
                println(io)
            end
        end
    end

    println("Done: $OUT_GRDECL")
    println("Z position fix: Lower Malm is the base; modelled deposits are exported above it.")
    println("Vertical shift applied upward: $(VERTICAL_SHIFT_UP_M) m")
    println("I/J handedness fix used: FLIP_X = false, FLIP_Y = true.")
    println("If orientation is mirrored the wrong way, try FLIP_X = true, FLIP_Y = false instead.")
    println("If the grid imports upside down vertically, set Z_MODE = :depth and rerun.")
end

export_grdecl()
