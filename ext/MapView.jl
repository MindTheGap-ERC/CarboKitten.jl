# ~/~ begin <<docs/src/visualization/map-view.md#ext/MapView.jl>>[init]
module MapView

import CarboKitten.Visualization: map_view, map_view!
using CarboKitten.Utility: in_units_of
using CarboKitten.Export: Header, Data, DataVolume, read_volume
using CarboKitten.Output.Abstract: stratigraphic_column, water_depth

using Makie
using Unitful

const Time = typeof(1.0u"Myr")

# Pulled verbatim from WheelerDiagram.dominant_facies! — the same trick works on
# any array whose first axis is facies. For a DataVolume snapshot the input is
# (n_facies, n_x, n_y) and the output is (n_x, n_y) Int.
_colormax(d::AbstractArray) = getindex.(argmax(d; dims=1)[1, :, :], 1)

# Resolve a stratigraphic position into an index along the (write-interval
# corrected) time axis. Integer indices are passed through; Unitful time
# quantities are matched by nearest neighbour.
_to_time_index(t_axis::AbstractVector, idx::Integer) = idx
_to_time_index(t_axis::AbstractVector{<:Quantity}, t::Quantity) =
    argmin(abs.(t_axis .- t))

"""
    map_view!(ax, header, data;
              time             = end,
              show             = :preserved,
              colors           = Makie.wong_colors(),
              mask_emerged     = true,
              show_shoreline   = false,
              shoreline_kwargs = (color = :black, linewidth = 1.5),
              kwargs...)

Plot a single map view of the platform onto `ax`, colored by dominant facies at
the chosen stratigraphic position.

# Arguments
- `ax::Makie.Axis` — target axis.
- `header::Header` — simulation metadata.
- `data::DataVolume` — volume output (3-D dataset in x, y, t).

# Keyword arguments
- `time` — stratigraphic position. Either an integer write-frame index (1-based;
  defaults to the final frame) or a `Unitful.Quantity` time value such as
  `0.5u"Myr"`, in which case the nearest available frame is used.
- `show::Symbol` — `:model` shows what is being deposited in the chosen frame;
  `:preserved` (default) shows what is preserved in the stratigraphic column for
  that frame; `:both` overlays a translucent `:model` layer underneath the
  `:preserved` one. Same semantics as `WheelerDiagram.dominant_facies!`.
- `colors` — vector of facies colors, defaults to `Makie.wong_colors()`.
- `mask_emerged::Bool` — if `true` (default), cells that are emerged at the
  chosen frame (water_depth > 0 in CarboKitten's convention, i.e. elevation
  above sea level) are masked white.
- `show_shoreline::Bool` — if `true`, overlays a contour of the sea-level
  intersection (water_depth = 0) at the chosen frame.
- `shoreline_kwargs` — named tuple forwarded to `contour!` for the shoreline.
- `kwargs...` — forwarded to `heatmap!`.

# Returns
The `Heatmap` object (use it for attaching a `Colorbar`).
"""
function map_view!(ax::Makie.Axis, header::Header, data::DataVolume;
                   time::Union{Integer,Quantity} = length(header.axes.t[1:data.write_interval:end]),
                   show::Symbol = :preserved,
                   colors = Makie.wong_colors(),
                   mask_emerged::Bool = true,
                   show_shoreline::Bool = false,
                   shoreline_kwargs = (color = :black, linewidth = 1.5),
                   kwargs...)

    if show ∉ (:model, :preserved, :both)
        error("`show` must be one of :model, :preserved, :both — got $(show)")
    end

    prec = 1e-8u"m"  # below this we treat accumulation as zero (mirrors WheelerDiagram)
    n_facies = size(data.production, 1)
    wi = data.write_interval
    t_axis = header.axes.t[1:wi:end]
    t_idx = _to_time_index(t_axis, time)

    if t_idx < 1 || t_idx > length(t_axis)
        error("resolved time index $(t_idx) is out of range 1:$(length(t_axis))")
    end

    xkm = header.axes.x |> in_units_of(u"km")
    ykm = header.axes.y |> in_units_of(u"km")
    colormap = cgrad(colors[1:n_facies], n_facies, categorical = true)
    colorrange = (0.5, n_facies + 0.5)

    # Mask of emerged cells at this frame, computed once and reused for both
    # the model-mode mask and the optional shoreline contour.
    wd = water_depth(header, data)[:, :, t_idx]

    function model_facies()
        m = Matrix{Union{Missing,Int}}(_colormax(data.deposition[:, :, :, t_idx]))
        if mask_emerged
            m[wd .> 0u"m"] .= missing
        end
        return m
    end

    function preserved_facies()
        sc_full = stratigraphic_column(data)
        sc = sc_full[:, :, :, t_idx]
        m = Matrix{Union{Missing,Int}}(_colormax(sc))
        # blank out cells with effectively no preserved accumulation at this frame
        acc = dropdims(sum(sc; dims = 1); dims = 1)
        m[acc .< prec] .= missing
        return m
    end

    # Merge defaults with user kwargs so caller's keys cleanly override ours
    # (a bare splat would raise `duplicate keyword argument`).
    base   = (colormap = colormap, colorrange = colorrange, nan_color = :white)
    user   = (; kwargs...)
    hm_kw  = merge(base, user)

    hm = if show == :model
        heatmap!(ax, xkm, ykm, model_facies(); hm_kw...)
    elseif show == :preserved
        heatmap!(ax, xkm, ykm, preserved_facies(); hm_kw...)
    else  # :both — translucent model under preserved
        model_kw = merge((alpha = 0.3,), hm_kw)
        pres_kw  = merge(hm_kw, (nan_color = :transparent,))
        heatmap!(ax, xkm, ykm, model_facies();     model_kw...)
        heatmap!(ax, xkm, ykm, preserved_facies(); pres_kw...)
    end

    if show_shoreline
        contour!(ax, xkm, ykm, wd |> in_units_of(u"m");
            levels = [0.0], shoreline_kwargs...)
    end

    ax.xlabel = "x [km]"
    ax.ylabel = "y [km]"
    ax.aspect = DataAspect()
    t_myr = t_axis[t_idx] |> in_units_of(u"Myr")
    ax.title  = "t = $(round(t_myr; digits = 3)) Myr"

    return hm
end

"""
    map_view(header, data;
             times     = [end],
             layout    = :auto,
             colorbar  = true,
             size      = nothing,
             kwargs...) -> Makie.Figure

    map_view(filename, group; kwargs...) -> Makie.Figure

Build a figure with one map-view panel per stratigraphic position in `times`.

# Arguments
- `header::Header`, `data::DataVolume` — in-memory data (e.g. from
  `MemoryOutput`).
- *or* `filename::AbstractString`, `group::Union{Symbol,AbstractString}` —
  path to an HDF5 file and the name of the volume dataset inside it. This form
  calls `CarboKitten.Export.read_volume(filename, group)` internally.

# Keyword arguments
- `times` — vector of stratigraphic positions. Each element may be an integer
  write-frame index or a `Unitful.Quantity` time value. Defaults to a single
  panel at the final frame.
- `layout::Symbol` — `:row`, `:col`, or `:auto` (default; chooses a near-square
  grid).
- `colorbar::Bool` — add a shared colorbar for facies (default `true`).
- `size` — Makie figure size tuple; if `nothing`, a size is chosen from the
  layout.
- All other `kwargs` are forwarded to `map_view!`. Notably: `show`,
  `mask_emerged`, `show_shoreline`, `colors`.

# Examples

```julia
using GLMakie, Unitful
using CarboKitten.Export: read_volume
using CarboKitten.Visualization: map_view, map_view!

# 1. From an HDF5 file, one panel per chosen time
fig = map_view("data/output/alcap-example.h5", :topography;
               times = [0.2u"Myr", 0.5u"Myr", 1.0u"Myr"],
               show = :preserved,
               show_shoreline = true)

# 2. From in-memory data (e.g. a MemoryOutput result)
fig = map_view(result.header, result.data_volumes[:topography];
               times = [10, 50, 100])
# ~/~ end
