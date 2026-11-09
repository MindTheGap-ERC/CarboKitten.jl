# ~/~ begin <<docs/src/visualization/fence-diagrams.md#ext/FenceDiagrams.jl>>[init]
module FenceDiagram

import CarboKitten.Visualization: fence_diagram, fence_diagram!

using CarboKitten.Visualization
using CarboKitten.Utility: in_units_of
using CarboKitten.Output.Abstract: AbstractBundle, load, load_volume, load_slice
using CarboKitten.Output.Storage: Header, DataSlice, DataVolume, Pack, PackSlice, PackVolume,
    header, data_slices, deposition, water_depth, surface_heights, n_facies, slice,
    time_axis, x_axis, y_axis, sea_level, bedrock
using CarboKitten.Output.MemoryWriter: MemoryOutput
using CarboKitten.Algorithms: skeleton

import ..SedimentProfile: explode_quad_vertices

using Makie
using GeometryBasics
using Unitful

const Time   = typeof(1.0u"Myr")

_to_index(axis::AbstractVector, idx::Integer) = Int(idx)
_to_index(axis::AbstractVector{<:Quantity}, pos::Quantity) =
    argmin(abs.(axis .- pos))

function _facies_fraction(column, facies::Integer)
    total = sum(column)
    iszero(total) && return NaN
    return column[facies] / total
end

function _slice_geometry(p::PackSlice)
    s = slice(p)
    if s[1] isa Colon && s[2] isa Integer
        return :along_x, x_axis(p), y_axis(p)
    elseif s[1] isa Integer && s[2] isa Colon
        return :along_y, y_axis(p), x_axis(p)
    else
        error("fence_diagram: unsupported slice $(s); need (:, Int) or (Int, :)")
    end
end

function fence_plot!(ax::Axis3, p::PackSlice;
                     color::AbstractArray, mesh_args...)
    orient, pos_axis, fixed_coord = _slice_geometry(p)

    pos_km   = pos_axis    |> in_units_of(u"km")
    fixed_km = fixed_coord |> in_units_of(u"km")
    h_m      = surface_heights(p) |> in_units_of(u"m")

    verts = zeros(Float64, size(h_m)..., 3)
    if orient === :along_x
        @views verts[:, :, 1] .= pos_km
        @views verts[:, :, 2] .= fixed_km
        @views verts[:, :, 3] .= h_m
    else
        @views verts[:, :, 1] .= fixed_km
        @views verts[:, :, 2] .= pos_km
        @views verts[:, :, 3] .= h_m
    end

    v, f = explode_quad_vertices(verts)
    c = vec(color)
    return mesh!(ax, v, f; color=vcat(c, c), mesh_args...)
end

function fence_plot!(f::F, ax::Axis3, p::PackSlice; mesh_args...) where {F}
    color = f.(eachslice(deposition(p), dims=(2, 3)))
    return fence_plot!(ax, p; color=color, mesh_args...)
end

# Bedrock and sea-level only need the Header — initial_topography, subsidence_rate,
# sea_level and axes are all stored there. No DataVolume required.

function _plot_bedrock!(ax::Axis3, hdr::Header; alpha::Real=0.35)
    x_km = x_axis(hdr) |> in_units_of(u"km")
    y_km = y_axis(hdr) |> in_units_of(u"km")
    base = bedrock(hdr) |> in_units_of(u"m")
    surface!(ax, x_km, y_km, base;
             color=fill(0.5, size(base)), colormap=:grays,
             alpha=alpha, transparency=true, colorrange=(0.0, 1.0))
end

function _plot_sealevel!(ax::Axis3, hdr::Header; alpha::Real=0.2)
    x_km = x_axis(hdr) |> in_units_of(u"km")
    y_km = y_axis(hdr) |> in_units_of(u"km")
    sl_end = last(sea_level(hdr)) |> in_units_of(u"m")
    z = fill(sl_end, length(x_km), length(y_km))
    surface!(ax, x_km, y_km, z;
             color=fill(0.7, size(z)), colormap=:Blues,
             alpha=alpha, transparency=true, colorrange=(0.0, 1.0))
end

function _plot_fence_unconformities!(ax::Axis3, p::PackSlice,
                                     h_m::AbstractMatrix; minwidth::Int, kwargs...)
    orient, pos_axis, fixed_coord = _slice_geometry(p)
    pos_km   = pos_axis    |> in_units_of(u"km")
    fixed_km = fixed_coord |> in_units_of(u"km")
    hiatus = skeleton(water_depth(p) .< 0.0u"m"; minwidth=minwidth)
    isempty(hiatus[1]) && return
    verts = if orient === :along_x
        [Point3f(pos_km[pt[1]], fixed_km, h_m[pt...]) for pt in hiatus[1]]
    else
        [Point3f(fixed_km, pos_km[pt[1]], h_m[pt...]) for pt in hiatus[1]]
    end
    linesegments!(ax, vec(permutedims(verts[hiatus[2]])); kwargs...)
end

function _plot_fence_coeval_lines!(ax::Axis3, p::PackSlice,
                                   h_m::AbstractMatrix, tics::Vector{Int}; kwargs...)
    orient, pos_axis, fixed_coord = _slice_geometry(p)
    pos_km   = pos_axis    |> in_units_of(u"km")
    fixed_km = fixed_coord |> in_units_of(u"km")
    n_pos = size(h_m, 1)
    for t in tics
        t = clamp(t, 1, size(h_m, 2))
        if orient === :along_x
            lines!(ax, pos_km, fill(fixed_km, n_pos), h_m[:, t]; kwargs...)
        else
            lines!(ax, fill(fixed_km, n_pos), pos_km, h_m[:, t]; kwargs...)
        end
    end
end

_apply_unconformities!(::Axis3, ::PackSlice, _, ::Nothing) = nothing
_apply_unconformities!(ax::Axis3, p::PackSlice, hm, flag::Bool) =
    flag && _plot_fence_unconformities!(ax, p, hm;
        minwidth=10, color=:white, linestyle=:dash, linewidth=1)
_apply_unconformities!(ax::Axis3, p::PackSlice, hm, minwidth::Int) =
    _plot_fence_unconformities!(ax, p, hm;
        minwidth=minwidth, color=:white, linestyle=:dash, linewidth=1)

function _apply_coeval_lines!(ax::Axis3, p::PackSlice, h_m, n_tics::Tuple{Int,Int})
    n_steps  = length(time_axis(p)) - 1
    n_major, n_minor = n_tics
    major = collect(div.(n_steps:n_steps:n_steps*n_major, n_major) .+ 1)
    minor = filter(t -> !(t in major),
                   collect(div.(n_steps:n_steps:n_steps*n_minor, n_minor) .+ 1))
    _plot_fence_coeval_lines!(ax, p, h_m, minor;
        color=:black, linewidth=0.6, linestyle=:dot)
    _plot_fence_coeval_lines!(ax, p, h_m, major;
        color=:black, linewidth=0.8, linestyle=:solid)
end
_apply_coeval_lines!(ax::Axis3, p::PackSlice, hm, flag::Bool) =
    flag && _apply_coeval_lines!(ax, p, hm, (4, 8))
_apply_coeval_lines!(ax::Axis3, p::PackSlice, hm, tics::Vector{Int}) =
    _plot_fence_coeval_lines!(ax, p, hm, tics;
        color=:black, linewidth=0.8, linestyle=:solid)
function _apply_coeval_lines!(ax::Axis3, p::PackSlice, h_m, tics::Vector{<:Time})
    t_axis = time_axis(p)
    idx = Int[searchsortedfirst(t_axis, ti) for ti in tics]
    _plot_fence_coeval_lines!(ax, p, h_m, idx;
        color=:black, linewidth=0.8, linestyle=:solid)
end

"""
    fence_diagram!(ax::Axis3, slices; kwargs...)
    fence_diagram!(ax::Axis3, header::Header, slices; kwargs...)

Core fence-diagram renderer. `slices` is any iterable of `PackSlice` objects
(or of `DataSlice` objects, when a `header` is given); each slice is rendered
as a vertical panel at its actual 3D position (inferred from its slice
specification).

Bedrock and sea-level surfaces are produced from the header alone — they do
not require a volume — so this method produces a complete plot from just a
sequence of slices. This makes it possible to plot slices read from different
HDF5 files, extracted from a `MemoryOutput`, or selected arbitrarily from a
volume.

Keyword arguments:

- `show_unconformities::Union{Nothing,Bool,Int}=true`
- `show_coeval_lines::Union{Bool,Tuple{Int,Int},Vector{Int},Vector{Time}}=false`
- `show_bedrock::Bool=true`
- `show_sealevel::Bool=false`
- `color_by::Symbol=:facies` — `:facies` or `:facies_fraction`
- `facies::Union{Nothing,Integer}=nothing` — required when `color_by=:facies_fraction`
- `colormap`, `alpha` — forwarded to `mesh!`
"""
function fence_diagram!(ax::Axis3, slices;
                        show_unconformities::Union{Nothing,Bool,Int}=true,
                        show_coeval_lines::Union{Bool,Tuple{Int,Int},Vector{Int},Vector{<:Time}}=false,
                        show_bedrock::Bool=true,
                        show_sealevel::Bool=false,
                        color_by::Symbol=:facies,
                        facies::Union{Nothing,Integer}=nothing,
                        colormap=nothing,
                        alpha::Real=1.0,
                        mesh_args...)
    isempty(slices) && error("fence_diagram!: no slices given.")
    hdr = header(first(slices))
    show_bedrock  && _plot_bedrock!(ax, hdr)
    show_sealevel && _plot_sealevel!(ax, hdr)

    plot_ref = nothing

    color_function, cmap, colorrange = if color_by === :facies
        (argmax, nothing, nothing)
    elseif color_by === :facies_fraction
        facies === nothing &&
            error("fence_diagram!: `facies` must be specified when `color_by=:facies_fraction`.")
        facies_idx = Int(facies)
        (col -> _facies_fraction(col, facies_idx),
         colormap === nothing ? :viridis : colormap,
         (0.0, 1.0))
    else
        error("fence_diagram!: `color_by` must be :facies or :facies_fraction, got :$(color_by)")
    end

    for p in slices
        n_f = n_facies(p)
        resolved_cmap = if color_by === :facies
            colormap === nothing ?
                cgrad(Makie.wong_colors()[1:n_f], n_f, categorical=true) : colormap
        else
            cmap
        end
        resolved_range = color_by === :facies ? (1, n_f) : colorrange

        plt = fence_plot!(color_function, ax, p;
                          colormap=resolved_cmap, colorrange=resolved_range,
                          alpha=alpha, mesh_args...)
        plot_ref = something(plot_ref, plt)

        h_m = surface_heights(p) |> in_units_of(u"m")
        _apply_coeval_lines!(ax, p, h_m, show_coeval_lines)
        _apply_unconformities!(ax, p, h_m, show_unconformities)
    end

    ax.xlabel = "x [km]"
    ax.ylabel = "y [km]"
    ax.zlabel = "depth [m]"
    ax.title  = "fence diagram"
    return plot_ref
end

"""
    fence_diagram!(ax::Axis3, volume::PackVolume; x_slices=[], y_slices=[], kwargs...)
    fence_diagram!(ax::Axis3, header::Header, data::DataVolume; x_slices=[], y_slices=[], kwargs...)

Convenience method: generates slices from `x_slices` and `y_slices` and
delegates entirely to the sequence method. Slice positions can be given as
indices or as coordinates. All kwargs are forwarded.
"""
function fence_diagram!(ax::Axis3, p::PackVolume;
                        x_slices::AbstractVector=Int[],
                        y_slices::AbstractVector=Int[],
                        kwargs...)
    isempty(x_slices) && isempty(y_slices) &&
        error("fence_diagram!: specify at least one slice via `x_slices` or `y_slices`.")

    x_idx = Int[_to_index(x_axis(p), pos) for pos in x_slices]
    y_idx = Int[_to_index(y_axis(p), pos) for pos in y_slices]

    slices = [[p[i, :] for i in x_idx]; [p[:, j] for j in y_idx]]
    return fence_diagram!(ax, slices; kwargs...)
end

fence_diagram!(ax::Axis3, header::Header, data::DataVolume; kwargs...) =
    fence_diagram!(ax, Pack(header, data); kwargs...)
fence_diagram!(ax::Axis3, header::Header, slices; kwargs...) =
    fence_diagram!(ax, [Pack(header, s) for s in slices]; kwargs...)

"""
    fence_diagram(volume::PackVolume; x_slices=[], y_slices=[], kwargs...)
    fence_diagram(slices; kwargs...)
    fence_diagram(header::Header, data; kwargs...)

Create a figure with a fence diagram and a legend or colorbar. The `data`
is either a `PackVolume`, from which slices are taken, or an iterable of
`PackSlice`. See `fence_diagram!` for keyword arguments.
"""
function fence_diagram(data;
                       size::Tuple{Int,Int}=(1400, 800),
                       azimuth::Real=-π/3,
                       elevation::Real=π/8,
                       color_by::Symbol=:facies,
                       facies::Union{Nothing,Integer}=nothing,
                       colormap=nothing,
                       kwargs...)
    n_f = n_facies(data isa PackVolume ? data : first(data))

    fig = Figure(size=size)
    ax  = Axis3(fig[1, 1])
    ax.azimuth   = azimuth
    ax.elevation = elevation

    fence_diagram!(ax, data;
        color_by=color_by, facies=facies, colormap=colormap, kwargs...)

    if color_by === :facies
        colors   = Makie.wong_colors()[1:n_f]
        elements = [PolyElement(color=colors[i]) for i in 1:n_f]
        labels   = ["Facies $(i)" for i in 1:n_f]
        Legend(fig[1, 2], elements, labels, "Facies")
    elseif color_by === :facies_fraction
        Colorbar(fig[1, 2];
            colormap = colormap === nothing ? :viridis : colormap,
            limits   = (0, 1),
            label    = "Proportion of facies $(facies)")
    end
    return fig
end

fence_diagram(header::Header, data::DataVolume; kwargs...) =
    fence_diagram(Pack(header, data); kwargs...)
fence_diagram(header::Header, slices; kwargs...) =
    fence_diagram([Pack(header, s) for s in slices]; kwargs...)

"""
    fence_diagram(filename, group; x_slices=[], y_slices=[], kwargs...)

Plot slices of the volume data set `group` stored in `filename`.
"""
function fence_diagram(filename::AbstractString,
                       group::Union{Symbol,AbstractString}; kwargs...)
    load(filename) do bundle
        fence_diagram(load_volume(bundle, Symbol(group)); kwargs...)
    end
end

"""
    fence_diagram(bundle::AbstractBundle; kwargs...)
    fence_diagram(filename; kwargs...)
    fence_diagram(output::MemoryOutput; kwargs...)

Plot all slice data sets found in the given bundle, file or output as a fence diagram.
"""
function fence_diagram(bundle::AbstractBundle; kwargs...)
    labels = sort(data_slices(bundle))
    isempty(labels) && error("fence_diagram: no slice data sets found")
    return fence_diagram([load_slice(bundle, l) for l in labels]; kwargs...)
end

fence_diagram(filename::AbstractString; kwargs...) =
    load(bundle -> fence_diagram(bundle; kwargs...), filename)

fence_diagram(output::MemoryOutput; kwargs...) =
    fence_diagram(load(output); kwargs...)

end  # module
# ~/~ end
