# export_carbokitten_to_petrel_grdecl.jl
# Fixed version: dimension normalization + Julia 1.12 printf + fraction export fix
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
const OUT_GRDECL  = "data/output/eclepens_carbokitten_petrel.GRDECL"

# Use :bbox for a non-rotated Petrel grid covering the TSurf extent.
# Use :pca if your Petrel grid is rotated along the main TSurf trend.
const GRID_MODE = :bbox       # :bbox or :pca

# The uploaded TSurf says ZPOSITIVE Elevation and contains negative Z values.
# Use :elevation to preserve those Petrel elevations.
# If Petrel imports the grid upside down, switch to :depth.
const Z_MODE = :elevation     # :elevation or :depth

# HDF5 output frame 1 is normally zero thickness. Keep false for Petrel export.
const INCLUDE_ZERO_FIRST_FRAME = false

# Write continuous fraction properties as well as dominant integer properties.
# This makes the GRDECL larger but much more useful in Petrel.
const WRITE_FRACTIONS = true

# Cells thinner than this are inactive.
const MIN_CELL_THICKNESS_M = 0.05

# Extra vertical length added to COORD pillars below the deepest ZCORN.
const PILLAR_EXTRA_M = 50.0

# Invert the HDF5 Y index when mapping to Petrel coordinates.
# Toggle this if the Petrel map appears north-south flipped.
const FLIP_Y = false

# IDW interpolation settings for the TSurf top surface.
const IDW_K = 12
const IDW_POWER = 2.0

# Property labels.
const PROD_NAMES = [
    "ooids",
    "corals",
    "mud",
    "peloids",
    "bioclasts_intraclasts",
    "oncoids",
]

const CLASS_NAMES = [
    "FA5 tidal flat / restricted platform",
    "FA3 ooid shoal complex",
    "FA3.4/FA2.2 coral-microbialite buildups",
    "FA4.4 interior platform mudstone",
    "FA4 back-shoal peloid-oncoid",
    "FA2 mid-ramp oncoid/bioclastic",
    "FA1 lower offshore mudstone/wackestone",
    "FA2.4 upper-offshore bioclastic wackestone",
    "fallback",
]

const PROD_FRAC_KEYS = ["P_OOIDS", "P_CORAL", "P_MUD", "P_PELO", "P_BIOIN", "P_ONCO"]
const CLS_FRAC_KEYS  = ["C_FA5", "C_FA3", "C_REEF", "C_FA44", "C_FA4", "C_FA2", "C_FA1", "C_F24", "C_FALLB"]

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
            jj = FLIP_Y ? (ny + 2 - j) : j
            X[j, i] = xs[i]
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
            jj = FLIP_Y ? (ny + 2 - j) : j
            p = center .+ V * [us[i], vs[jj]]
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

function z_at(top_z::Float64, cumulative_thickness::Float64)
    if Z_MODE == :elevation
        return top_z - cumulative_thickness
    elseif Z_MODE == :depth
        return -top_z + cumulative_thickness
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
    nt, ny, nx = size(cumulative)
    nz = length(kframes)

    zcorn = zeros(Float64, 2*nx, 2*ny, 2*nz)

    prev_corners = centers_to_corners(zeros(Float64, ny, nx))

    for kk in 1:nz
        kt = kframes[kk]
        bottom_corners = centers_to_corners(cumulative[kt, :, :])

        for j in 1:ny, i in 1:nx
            top_sw = z_at(top_surface[j,   i],   prev_corners[j,   i])
            top_se = z_at(top_surface[j,   i+1], prev_corners[j,   i+1])
            top_nw = z_at(top_surface[j+1, i],   prev_corners[j+1, i])
            top_ne = z_at(top_surface[j+1, i+1], prev_corners[j+1, i+1])

            bot_sw = z_at(top_surface[j,   i],   bottom_corners[j,   i])
            bot_se = z_at(top_surface[j,   i+1], bottom_corners[j,   i+1])
            bot_nw = z_at(top_surface[j+1, i],   bottom_corners[j+1, i])
            bot_ne = z_at(top_surface[j+1, i+1], bottom_corners[j+1, i+1])

            zcorn[2*i-1, 2*j-1, 2*kk-1] = top_sw
            zcorn[2*i,   2*j-1, 2*kk-1] = top_se
            zcorn[2*i-1, 2*j,   2*kk-1] = top_nw
            zcorn[2*i,   2*j,   2*kk-1] = top_ne

            zcorn[2*i-1, 2*j-1, 2*kk] = bot_sw
            zcorn[2*i,   2*j-1, 2*kk] = bot_se
            zcorn[2*i-1, 2*j,   2*kk] = bot_nw
            zcorn[2*i,   2*j,   2*kk] = bot_ne
        end

        prev_corners = bottom_corners
    end

    return vec(zcorn)
end

function build_coord(X::Matrix{Float64}, Y::Matrix{Float64}, top_surface::Matrix{Float64},
                     cumulative::Array{Float64, 3})
    ny1, nx1 = size(X)
    nt, ny, nx = size(cumulative)

    max_cum_center = cumulative[end, :, :]
    max_cum_corner = centers_to_corners(max_cum_center)

    coord = Float64[]

    for j in 1:ny1, i in 1:nx1
        ztop = z_at(top_surface[j, i], 0.0)
        zbot = z_at(top_surface[j, i], max_cum_corner[j, i] + PILLAR_EXTRA_M)
        append!(coord, (X[j, i], Y[j, i], ztop, X[j, i], Y[j, i], zbot))
    end

    return coord
end

function build_actnum(cumulative::Array{Float64, 3}, kframes::Vector{Int})
    nt, ny, nx = size(cumulative)
    vals = Int[]

    prev = zeros(Float64, ny, nx)

    for kt in kframes
        cur = cumulative[kt, :, :]
        thick = cur .- prev
        for j in 1:ny, i in 1:nx
            push!(vals, thick[j, i] > MIN_CELL_THICKNESS_M ? 1 : 0)
        end
        prev = cur
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

    kframes = INCLUDE_ZERO_FIRST_FRAME ? collect(1:nt) : collect(2:nt)
    nz = length(kframes)

    println("Normalized internal sizes:")
    println("  production = $(size(prod))  interpreted as (nt, ny, nx, nprod)")
    println("  cumulative = $(size(cum))   interpreted as (nt, ny, nx)")
    println("  classified = $(size(cls))   interpreted as (nt, ny, nx, nclass)")
    println("Grid dimensions exported: nx=$nx ny=$ny nz=$nz")
    println("Production factories: $nprod")
    println("Classified classes: $ncls")

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

    println("Writing $OUT_GRDECL ...")
    mkpath(dirname(OUT_GRDECL))

    open(OUT_GRDECL, "w") do io
        println(io, "-- Exported from CarboKitten HDF5")
        println(io, "-- Production HDF5: $PROD_H5")
        println(io, "-- Classified HDF5: $CLASS_H5")
        println(io, "-- TSurf geometry: $TSURF_FILE")
        println(io, "-- GRID_MODE=$GRID_MODE Z_MODE=$Z_MODE")
        println(io, "--")
        println(io, "-- PRODFAC values:")
        for (i, name) in enumerate(PROD_NAMES)
            println(io, "--   $i = $name")
        end
        println(io, "--")
        println(io, "-- CLSFAC values:")
        for (i, name) in enumerate(CLASS_NAMES)
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

        println(io, "PRODFAC")
        write_int_values(io, prod_dom; perline=30)
        println(io)

        println(io, "CLSFAC")
        write_int_values(io, cls_dom; perline=30)
        println(io)

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
    println("If the grid is north-south flipped in Petrel, set FLIP_Y = true and rerun.")
    println("If the grid imports upside down vertically, set Z_MODE = :depth and rerun.")
end

export_grdecl()
