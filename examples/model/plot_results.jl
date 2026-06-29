# examples/model/plot_results.jl
#
# Reads the saved production and classified HDF5 files and plots:
#   - Fence diagrams production + classified
#   - Final map views production + classified
#   - Depth-slice map views production + classified
#   - Stratigraphic columns production + classified
#
# Usage:
#   julia --project=. examples/model/plot_results.jl

module PlotResults

using CarboKitten
using CarboKitten.Export: read_volume, read_slice
using CarboKitten.Visualization: fence_diagram!
using CarboKitten.Output.Abstract: stratigraphic_column
using GLMakie
using Unitful

include(joinpath(@__DIR__, "../../ext/MapView.jl"))
include(joinpath(@__DIR__, "../../ext/StratigraphicColumn.jl"))

using .MapView: map_view
using .StratigraphicColumn: stratigraphic_column!

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

const OUTPUT_FILE     = "data/output/eclepens_withoutca.h5"
const OUTPUT_CLS_FILE = "data/output/eclepens_withoutca_classified.h5"
const FIG_DIR         = "data/output/eclepens_withoutca_figs"

mkpath(FIG_DIR)

# ---------------------------------------------------------------------------
# Depths for depth-slice maps
# ---------------------------------------------------------------------------

const TARGET_DEPTHS_M = [50.0, 100.0, 150.0, 200.0]

# ---------------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------------

const PROD_LABELS = [
    "ooids",
    "corals",
    "mud",
    "peloids",
    "grains",
    "oncoids",
]

const CLS_LABELS = [
    "FA5 tidal flat / local exposure",
    "FA3 ooid shoal complex",
    "FA3/FA2 coral-microbialite buildup",
    "FA4 back-shoal peloid-oncoid",
    "FA4 interior platform mudstone",
    "FA1 lower offshore mudstone",
    "FA2 upper offshore / oncoid-bioclastic wackestone",
    "fallback",
]

# ---------------------------------------------------------------------------
# Palette
# ---------------------------------------------------------------------------

function palette(n)
    base = [
        RGBAf(0.000, 0.447, 0.698, 1.0),
        RGBAf(0.902, 0.624, 0.000, 1.0),
        RGBAf(0.000, 0.620, 0.451, 1.0),
        RGBAf(0.835, 0.369, 0.000, 1.0),
        RGBAf(0.800, 0.475, 0.655, 1.0),
        RGBAf(0.941, 0.894, 0.259, 1.0),
        RGBAf(0.337, 0.706, 0.914, 1.0),
        RGBAf(0.550, 0.337, 0.294, 1.0),
        RGBAf(0.600, 0.600, 0.600, 1.0),
        RGBAf(0.650, 0.800, 0.350, 1.0),
        RGBAf(0.950, 0.500, 0.500, 1.0),
        RGBAf(0.500, 0.300, 0.700, 1.0),
        RGBAf(0.300, 0.750, 0.750, 1.0),
        RGBAf(0.900, 0.700, 0.200, 1.0),
        RGBAf(0.400, 0.400, 0.800, 1.0),
        RGBAf(0.200, 0.500, 0.300, 1.0),
        RGBAf(0.100, 0.100, 0.100, 1.0),
    ]

    @assert n <= length(base) "Need at least $n colors"
    return base[1:n]
end

# ---------------------------------------------------------------------------
# Read production output
# ---------------------------------------------------------------------------

@info "Reading production output..."

header, vol  = read_volume(OUTPUT_FILE, :topography)
_,      prof = read_slice(OUTPUT_FILE,  :profile)

nx, ny = header.grid_size

prod_colors = palette(length(PROD_LABELS))

# ---------------------------------------------------------------------------
# Read classified output
# ---------------------------------------------------------------------------

@info "Reading classified output..."

cls_header, cls_vol = read_volume(OUTPUT_CLS_FILE, :classified)

n_cls = cls_header.n_facies

@assert n_cls <= length(CLS_LABELS) "Classified file has $n_cls facies but only $(length(CLS_LABELS)) labels were provided"

cls_colors = palette(n_cls)
cls_labels = CLS_LABELS[1:n_cls]

# Extract classified profile at same y-index used by the model profile
cls_prof = cls_vol[:, div(cls_header.grid_size[2], 2)]

# ---------------------------------------------------------------------------
# Build fence slices
# ---------------------------------------------------------------------------

slices = vcat(
    [vol[i, :] for i in [div(nx, 4), div(nx, 2), 3 * div(nx, 4)] if i >= 1],
    [vol[:, j] for j in [div(ny, 4), div(ny, 2), 3 * div(ny, 4)] if j >= 1],
)

cls_slices = vcat(
    [cls_vol[i, :] for i in [div(nx, 4), div(nx, 2), 3 * div(nx, 4)] if i >= 1],
    [cls_vol[:, j] for j in [div(ny, 4), div(ny, 2), 3 * div(ny, 4)] if j >= 1],
)

@info "$(length(slices)) production slices"
@info "$(length(cls_slices)) classified slices"

# Stratigraphic column locations along profile
n_prof       = size(prof.sediment_thickness, 1)
sc_locs      = [div(n_prof, 4), div(n_prof, 2), 3 * div(n_prof, 4)]
phys_scale_m = ustrip(u"m", header.axes.x[2] - header.axes.x[1])

# ---------------------------------------------------------------------------
# Fence diagrams
# ---------------------------------------------------------------------------

function named_fence(header, slices, colors, labels; kwargs...)
    n_f = length(labels)

    fig = Figure(size=(1400, 650))
    ax  = Axis3(fig[1, 1])

    fence_diagram!(
        ax,
        header,
        slices;
        colormap = cgrad(colors[1:n_f], n_f; categorical=true),
        kwargs...,
    )

    Legend(
        fig[1, 2],
        [PolyElement(color=colors[i]) for i in 1:n_f],
        labels;
        framevisible=false,
    )

    return fig
end

@info "Plotting production fence diagram..."

fig_fence_prod = named_fence(
    header,
    slices,
    prod_colors,
    PROD_LABELS;
    color_by            = :facies,
    show_unconformities = 10,
    show_coeval_lines   = (1, 5),
    show_sealevel       = true,
)

save(joinpath(FIG_DIR, "fence_production.png"), fig_fence_prod; px_per_unit=2)

@info "Plotting classified fence diagram..."

fig_fence_cls = named_fence(
    cls_header,
    cls_slices,
    cls_colors,
    cls_labels;
    color_by            = :facies,
    show_unconformities = 10,
    show_coeval_lines   = (1, 5),
    show_sealevel       = true,
)

save(joinpath(FIG_DIR, "fence_classified.png"), fig_fence_cls; px_per_unit=2)

# ---------------------------------------------------------------------------
# Final map views
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Final map views
# ---------------------------------------------------------------------------

function legend_banks(n_labels)
    return n_labels <= 6 ? 2 : 3
end

function named_map(header, vol, colors, labels, frame_id; kwargs...)
    n_f = length(labels)

    fig = map_view(
        header,
        vol;
        times    = [frame_id],
        colorbar = false,
        colormap = cgrad(colors[1:n_f], n_f; categorical=true),
        kwargs...,
    )

    # Make the map figure larger. This prevents long classified labels from
    # squeezing the map out of the visible window.
    resize!(fig.scene, (1500, 950))

    # Put legend below the map, not on the right.
    Legend(
        fig[2, 1],
        [PolyElement(color=colors[i]) for i in 1:n_f],
        labels;
        orientation  = :horizontal,
        nbanks       = legend_banks(n_f),
        framevisible = false,
        tellwidth    = false,
    )

    rowgap!(fig.layout, 8)

    return fig
end

n_frames     = length(header.axes.t[1:vol.write_interval:end])
n_cls_frames = length(cls_header.axes.t[1:cls_vol.write_interval:end])

@info "Plotting production final map view..."

fig_map_prod = named_map(
    header,
    vol,
    prod_colors,
    PROD_LABELS,
    n_frames;
    color_by       = :facies,
    show           = :preserved,
    show_shoreline = false,
    mask_emerged   = true,
    layout         = :row,
)

save(joinpath(FIG_DIR, "mapview_production.png"), fig_map_prod; px_per_unit=2)

@info "Plotting classified final map view..."

fig_map_cls = named_map(
    cls_header,
    cls_vol,
    cls_colors,
    cls_labels,
    n_cls_frames;
    color_by       = :facies,
    show           = :preserved,
    show_shoreline = false,
    mask_emerged   = true,
    layout         = :row,
)

save(joinpath(FIG_DIR, "mapview_classified.png"), fig_map_cls; px_per_unit=2)

# ---------------------------------------------------------------------------
# Depth-slice map views
#
# For each cell, this finds the facies occupying a given depth below the final
# preserved surface. If the preserved column is thinner than the requested
# depth, the cell is plotted as light grey.
# ---------------------------------------------------------------------------

@info "Plotting depth-slice map views..."

function facies_at_depth(vol, depth_below_surface_m)
    # stratigraphic_column(vol):
    #   array dimensions are approximately:
    #   n_facies, nx, ny, n_t
    #   preserved net deposition, oldest to youngest
    sc = stratigraphic_column(vol)

    n_facies, nx, ny, n_t = size(sc)
    target = depth_below_surface_m * u"m"

    result = fill(NaN, nx, ny)

    for ix in 1:nx, iy in 1:ny
        layers = dropdims(sum(sc[:, ix, iy, :]; dims=1); dims=1)

        depth = zero(target)

        for t in n_t:-1:1
            depth += layers[t]

            if depth >= target
                dep = sc[:, ix, iy, t]
                f   = argmax(dep)

                if ustrip(dep[f]) > 0
                    result[ix, iy] = Float64(f)
                end

                break
            end
        end
    end

    return result
end

function plot_facies_at_depth(header, vol, colors, labels, depth_m, filename; title_prefix="Facies")
    n_f  = length(labels)
    xkm  = ustrip.(u"km", header.axes.x)
    ykm  = ustrip.(u"km", header.axes.y)
    cmap = cgrad(colors[1:n_f], n_f; categorical=true)

    # Larger window, with map on top and legend below.
    fig = Figure(size=(1500, 950))

    ax = Axis(
        fig[1, 1];
        title  = "$(title_prefix) at $(round(Int, depth_m)) m below final surface",
        xlabel = "x [km]",
        ylabel = "y [km]",
    )

    depth_slice = facies_at_depth(vol, depth_m)

    heatmap!(
        ax,
        xkm,
        ykm,
        depth_slice;
        colormap   = cmap,
        colorrange = (0.5, n_f + 0.5),
        nan_color  = :lightgrey,
    )

    ax.aspect = DataAspect()

    # Force the full model window to be visible.
    xlims!(ax, minimum(xkm), maximum(xkm))
    ylims!(ax, minimum(ykm), maximum(ykm))

    Legend(
        fig[2, 1],
        [PolyElement(color=colors[i]) for i in 1:n_f],
        labels;
        orientation  = :horizontal,
        nbanks       = legend_banks(n_f),
        framevisible = false,
        tellwidth    = false,
    )

    rowsize!(fig.layout, 1, Relative(0.82))
    rowsize!(fig.layout, 2, Relative(0.18))
    rowgap!(fig.layout, 8)

    save(filename, fig; px_per_unit=2)
    @info "  -> $(filename)"

    return fig
end

for depth_m in TARGET_DEPTHS_M
    @info "Depth slice at $(depth_m) m"

    plot_facies_at_depth(
        header,
        vol,
        prod_colors,
        PROD_LABELS,
        depth_m,
        joinpath(FIG_DIR, "mapview_depth$(round(Int, depth_m))m_production.png");
        title_prefix = "Production facies",
    )

    plot_facies_at_depth(
        cls_header,
        cls_vol,
        cls_colors,
        cls_labels,
        depth_m,
        joinpath(FIG_DIR, "mapview_depth$(round(Int, depth_m))m_classified.png");
        title_prefix = "Classified environment",
    )
end

# ---------------------------------------------------------------------------
# Stratigraphic columns
# ---------------------------------------------------------------------------

@info "Plotting production stratigraphic columns..."

fig_sc_prod = Figure(size=(300 * length(sc_locs) + 120, 600))

for (k, loc) in enumerate(sc_locs)
    ax = Axis(
        fig_sc_prod[1, k];
        title  = "x=$(round(Int, loc * phys_scale_m))m",
        xlabel = "depth [m]",
    )

    stratigraphic_column!(ax, header, prof[loc]; color=prod_colors)

    k > 1 && hideydecorations!(ax; ticks=false)
end

Legend(
    fig_sc_prod[1, end + 1],
    [PolyElement(color=prod_colors[i]) for i in 1:length(PROD_LABELS)],
    PROD_LABELS,
    "Factory";
    framevisible=false,
)

save(joinpath(FIG_DIR, "strat_columns_production.png"), fig_sc_prod; px_per_unit=2)

@info "Plotting classified stratigraphic columns..."

fig_sc_cls = Figure(size=(300 * length(sc_locs) + 150, 600))

for (k, loc) in enumerate(sc_locs)
    ax = Axis(
        fig_sc_cls[1, k];
        title  = "x=$(round(Int, loc * phys_scale_m))m",
        xlabel = "depth [m]",
    )

    stratigraphic_column!(ax, cls_header, cls_prof[loc]; color=cls_colors)

    k > 1 && hideydecorations!(ax; ticks=false)
end

Legend(
    fig_sc_cls[1, end + 1],
    [PolyElement(color=cls_colors[i]) for i in 1:n_cls],
    cls_labels,
    "Environment";
    framevisible=false,
)

save(joinpath(FIG_DIR, "strat_columns_classified.png"), fig_sc_cls; px_per_unit=2)

@info "Done. All figures in $(FIG_DIR)/"

end  # module PlotResults
