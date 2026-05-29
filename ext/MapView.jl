# ~/~ begin <<docs/src/visualization/map-view.md#ext/MapView.jl>>[init]
module MapView

import CarboKitten.Visualization: map_view, map_view!
using CarboKitten.Utility: in_units_of
using CarboKitten.Export: Header, Data, DataVolume, read_volume
using CarboKitten.Output.Abstract: stratigraphic_column, water_depth

using Makie
using Unitful

const Time = typeof(1.0u"Myr")

# Pulled verbatim from WheelerDiagram.dominant_facies! — works on
# any array whose first axis is facies. For a DataVolume snapshot the input is
# (n_facies, n_x, n_y) and the output is (n_x, n_y) Int.
_colormax(d::AbstractArray) = getindex.(argmax(d; dims=1)[1, :, :], 1)

#Calculate a given facies' proportion relative to the others in a specific location.
function _facies_fraction(d::AbstractArray, facies::Integer)
    total = dropdims(sum(d; dims = 1); dims = 1)
    selected = d[facies, :, :]

    fraction = Matrix{Union{Missing,Float64}}(undef, Base.size(total))

    for I in eachindex(total)
        if iszero(total[I])
            fraction[I] = missing
        else
            fraction[I] = Float64(ustrip(selected[I] / total[I]))
        end
    end

    return fraction
end

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
              color_by         = :facies,
              facies           = nothing,
              colors           = Makie.wong_colors(),
              colormap         = nothing,
              mask_emerged     = true,
              show_shoreline   = false,
              shoreline_kwargs = (color = :black, linewidth = 1.5),
              kwargs...)

Plot a single map view of the platform onto `ax`.

The map can be coloured either by dominant facies using categorical colours, or
by the proportion of a selected facies using a continuous colour scale..

# Arguments
- `ax::Makie.Axis` — target axis.
- `header::Header` — simulation metadata.
- `data::DataVolume` — volume output (3-D dataset in x, y, t).

# Keyword arguments
- `time` — stratigraphic position. Either an integer write-frame index
  (1-based; defaults to the final frame) or a `Unitful.Quantity` time value
  such as `0.5u"Myr"`, in which case the nearest available frame is used.
- `show::Symbol` — controls which stratigraphic quantity is plotted.
  `:model` shows what is being deposited in the chosen frame; `:preserved`
  shows what is preserved in the stratigraphic column for that frame; `:both`
  overlays a translucent `:model` layer underneath the `:preserved` layer.
- `color_by::Symbol` — controls how cells are coloured. Use `:facies` for
  categorical dominant-facies colouring, or `:facies_fraction` to colour by
  the proportion of one selected facies.
- `facies::Union{Nothing,Integer}` — facies index used when
  `color_by = :facies_fraction`.
- `colors` — vector of categorical facies colours used when
  `color_by = :facies`. Defaults to `Makie.wong_colors()`.
- `colormap` — colormap used for plotting. If `color_by = :facies`, this
  overrides the categorical facies colormap. If `color_by = :facies_fraction`,
  this controls the continuous colour scale and defaults to `:viridis`.
- `mask_emerged::Bool` — if `true`, cells that are emerged at the chosen frame
  are masked white.
- `show_shoreline::Bool` — if `true`, overlays a contour of the sea-level
  intersection (`water_depth = 0`) at the chosen frame.
- `shoreline_kwargs` — named tuple forwarded to `contour!` for the shoreline.
- `kwargs...` — forwarded to `heatmap!`.

# Returns
The `Heatmap` object. This can be used to attach a `Colorbar`.
"""

function map_view!(ax::Makie.Axis, header::Header, data::DataVolume;
                   time::Union{Integer,Quantity} = length(header.axes.t[1:data.write_interval:end]),
                   show::Symbol = :preserved,
                   color_by::Symbol = :facies,
                   facies::Union{Nothing,Integer} = nothing,
                   colors = Makie.wong_colors(),
                   colormap = nothing,
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

    if color_by == :facies
        cmap = colormap === nothing ?
            cgrad(colors[1:n_facies], n_facies, categorical = true) :
            colormap
        colorrange = (0.5, n_facies + 0.5)
    
    elseif color_by == :facies_fraction
        facies === nothing &&
            error("map_view!: `facies` must be specified when `color_by = :facies_fraction`.")
    
        facies_idx = Int(facies)
    
        if facies_idx < 1 || facies_idx > n_facies
            error("map_view!: `facies` must be between 1 and $(n_facies).")
        end
    
        cmap = colormap === nothing ? :viridis : colormap
        colorrange = (0.0, 1.0)
    
    else
        error("map_view!: `color_by` must be either :facies or :facies_fraction.")
    end

    # Mask of emerged cells at this frame, computed once and reused for both
    # the model-mode mask and the optional shoreline contour.
    wd = water_depth(header, data)[:, :, t_idx]

    function model_values()
        d = data.deposition[:, :, :, t_idx]
    
        m = if color_by == :facies
            Matrix{Union{Missing,Int}}(_colormax(d))
        else
            _facies_fraction(d, facies_idx)
        end
    
        if mask_emerged
            m[wd .> 0u"m"] .= missing
        end
    
        return m
    end
    
    function preserved_values()
        sc_full = stratigraphic_column(data)
        sc = sc_full[:, :, :, t_idx]
    
        m = if color_by == :facies
            Matrix{Union{Missing,Int}}(_colormax(sc))
        else
            _facies_fraction(sc, facies_idx)
        end
    
        acc = dropdims(sum(sc; dims = 1); dims = 1)
        m[acc .< prec] .= missing
    
        return m
    end

    # Merge defaults with user kwargs so caller's keys cleanly override ours
    base   = (colormap = cmap, colorrange = colorrange, nan_color = :white)
    user   = (; kwargs...)
    hm_kw  = merge(base, user)

    hm = if show == :model
        heatmap!(ax, xkm, ykm, model_values(); hm_kw...)
    elseif show == :preserved
        heatmap!(ax, xkm, ykm, preserved_values(); hm_kw...)
    else  # :both — translucent model under preserved
        model_kw = merge((alpha = 0.3,), hm_kw)
        pres_kw  = merge(hm_kw, (nan_color = :transparent,))
        heatmap!(ax, xkm, ykm, model_values();     model_kw...)
        heatmap!(ax, xkm, ykm, preserved_values(); pres_kw...)
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
             color_by::Symbol = :facies,
             facies::Union{Nothing,Integer} = nothing,
             colormap = nothing,
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
- `layout::Symbol` — `:row`, `:col`, or `:auto`. The default chooses a
  near-square grid.
- `colorbar::Bool` — add a shared colourbar or legend-like colour scale
  for the plotted maps. For `color_by = :facies`, the colourbar is categorical
  and labelled by facies index. For `color_by = :facies_fraction`, the colourbar
  is continuous from 0 to 1.
- `size` — Makie figure size tuple. If `nothing`, a size is chosen from the
  layout.
- `color_by::Symbol` — controls how cells are coloured. Use `:facies` for
  categorical dominant-facies colouring, or `:facies_fraction` to colour by
  the proportion of one selected facies.
- `facies::Union{Nothing,Integer}` — facies index used when
  `color_by = :facies_fraction`.
- `colormap` — optional colormap override. Defaults to categorical Wong colours
  for `color_by = :facies` and to `:viridis` for `color_by = :facies_fraction`.
- All other `kwargs` are forwarded to `map_view!`. Notably: `show`,
  `mask_emerged`, `show_shoreline`, `shoreline_kwargs`, and `colors`.
"""
function map_view(header::Header, data::DataVolume;
                  times::AbstractVector = [length(header.axes.t[1:data.write_interval:end])],
                  layout::Symbol = :auto,
                  colorbar::Bool = true,
                  size = nothing,
                  color_by::Symbol = :facies,
                  facies::Union{Nothing,Integer} = nothing,
                  colormap = nothing,
                  kwargs...)

    n = length(times)
    if n == 0
        error("`times` must contain at least one stratigraphic position")
    end

    nrows, ncols = if layout == :row
        (1, n)
    elseif layout == :col
        (n, 1)
    elseif layout == :auto
        ncols = ceil(Int, sqrt(n))
        nrows = ceil(Int, n / ncols)
        (nrows, ncols)
    else
        error("`layout` must be :row, :col, or :auto — got $(layout)")
    end

    fig_size = something(size, (320 * ncols + (colorbar ? 120 : 40), 320 * nrows + 60))
    fig = Figure(size = fig_size)

    hm_ref = nothing
    for (i, t) in enumerate(times)
        r = div(i - 1, ncols) + 1
        c = mod(i - 1, ncols) + 1
        ax = Axis(fig[r, c])
        hm = map_view!(ax, header, data;
            time = t,
            color_by = color_by,
            facies = facies,
            colormap = colormap,
            kwargs...)
        hm_ref = something(hm_ref, hm)
    end

    if colorbar && hm_ref !== nothing
        if color_by == :facies
            n_facies = Base.size(data.production, 1)
            Colorbar(fig[1:nrows, ncols + 1], hm_ref;
                     ticks = 1:n_facies,
                     label = "dominant facies")
        elseif color_by == :facies_fraction
            Colorbar(fig[1:nrows, ncols + 1], hm_ref;
                     ticks = 0:0.25:1,
                     label = "proportion of facies $(facies)")
        end
    end

    return fig
end

function map_view(filename::AbstractString,
                  group::Union{Symbol,AbstractString};
                  kwargs...)
    header, data = read_volume(filename, group)
    return map_view(header, data; kwargs...)
end

end  # module MapView
# ~/~ end
