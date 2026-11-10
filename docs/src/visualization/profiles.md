# Profile Plotting

![Sediment profile](../fig/sediment_profile.png)

The sediment profile is probably the most important visualization that we provide. By default it allows us to study the sediment composition of a section, by plotting the `argmax` of the deposition. In cases where significant amounts of sediment is eroded, all deposition is plotted, and it is assumed that newest depositions are shown on top of possible older ones.

``` {.julia .task file=examples/visualization/sediment_profile.jl}
#| creates: docs/src/_fig/sediment_profile.png
#| requires: data/output/alcap-example.h5
#| collect: figures

module Script
using CairoMakie
using CarboKitten.Export: read_slice
using CarboKitten.Visualization: sediment_profile

function main()
    save("docs/src/_fig/sediment_profile.png",
        sediment_profile(read_slice("data/output/alcap-example.h5", :profile)...))
end
end

Script.main()
```

If you want to visualize something other than the `argmax` of the deposition, you may use the `profile_plot!` function. For example, we can plot the fraction of second facies over total.

![](../fig/profile_fraction.png)

``` {.julia .task file=examples/visualization/profile_fraction.jl}
#| creates: docs/src/_fig/profile_fraction.png
#| requires: data/output/alcap-example.h5
#| collect: figures

module Script
using GLMakie
using CarboKitten.Export: read_slice
using CarboKitten.Visualization: profile_plot!

function main()
    (header, slice) = read_slice("data/output/alcap-example.h5", :profile)
    fig = Figure()
    ax = Axis(fig[1, 1])

    x = header.axes.x
    t = header.axes.t

    plot = profile_plot!(x -> x[2]/sum(x), ax, header, slice; colorrange=(0, 1))
    Colorbar(fig[1, 2], plot; label=L"f_2 / f_{total}")

    save("docs/src/_fig/profile_fraction.png", fig)
    fig
end
end

Script.main()
```

## Cross-section direction

A profile `DataSlice` can represent either a **dip section** (varying along x, fixed y index) or a **strike section** (fixed x index, varying along y). Both are produced by the same `DataVolume` indexing:

```julia
vol[:, j]   # dip section at y = j
vol[i, :]   # strike section at x = i
```

The correct spatial positions are inferred automatically from the slice specification (see `section_axis`) — no direction argument is needed. If `slice(pack)[1]` is an `Int`, the section is fixed in x and varies along y (strike); otherwise it varies along x (dip).

``` {.julia #section-positions}
"""
    _section_positions(pack)

Return the physical positions (in km) along the spatial axis of a profile.
For a dip section (varying x, fixed y) this is the x axis; for a strike
section (fixed x, varying y) this is the y axis. See `section_axis`.
"""
_section_positions(p::PackSlice) = section_axis(p) |> in_units_of(u"km")
```

## Proportion plot

Beyond the dominant-facies categorical view, it is useful to see the proportion of a specific facies relative to total sediment at each location and time step. The `sediment_proportion!` function provides this using the same mesh geometry as `sediment_profile!`.

![Proportion plot](../fig/profile_fraction.png)

``` {.julia .task file=examples/visualization/profile_proportion.jl}
#| creates: docs/src/_fig/profile_proportion.png
#| requires: data/output/alcap-example.h5
#| collect: figures

module Script
using CairoMakie
using CarboKitten.Export: read_slice
using CarboKitten.Visualization: sediment_proportion

function main()
    header, data = read_slice("data/output/alcap-example.h5", :profile)
    save("docs/src/_fig/profile_proportion.png",
         sediment_proportion(header, data, 1))
end
end

Script.main()
```

## Cross-section examples

The following script uses the existing `alcap-example.h5` volume output (the
`:topography` group, which carries full 2D spatial data at every 10 write steps).
It takes one dip section and one strike section from the mid-grid position, then
plots all four combinations of section direction × proportion mode.

``` {.julia .task file=examples/visualization/cross_sections.jl}
#| creates: docs/src/_fig/xsec_dip.png
#|          docs/src/_fig/xsec_strike.png
#|          docs/src/_fig/xsec_proportion_deposited.png
#|          docs/src/_fig/xsec_proportion_preserved.png
#| requires: data/output/alcap-example.h5
#| collect: figures

module Script

using CairoMakie
using CarboKitten.Export: read_volume
using CarboKitten.Visualization: sediment_profile, sediment_proportion

function main()
    header, vol = read_volume("data/output/alcap-example.h5", :topography)

    nx, ny = size(vol.sediment_thickness)[1:2]
    dip    = vol[:, div(ny, 2) + 1]   # dip section at mid-y
    strike = vol[div(nx, 2) + 1, :]   # strike section at mid-x

    save("docs/src/_fig/xsec_dip.png",
         sediment_profile(header, dip))

    save("docs/src/_fig/xsec_strike.png",
         sediment_profile(header, strike))

    save("docs/src/_fig/xsec_proportion.png",
         sediment_proportion(header, dip, 1))
end

end

Script.main()
```

Dip section (perpendicular to strike), dominant deposited facies:

![Dip cross section](../fig/xsec_dip.png)

Strike section (parallel to platform margin), dominant deposited facies:

![Strike cross section](../fig/xsec_strike.png)

Proportion of facies 1:

![Proportion deposited](../fig/xsec_proportion.png)

### Exploding Vertices

The profile plot is rendered using Makie's `mesh` function. This function accepts an array of (two-dimensional) vertices and an array of faces, usually triangles. The way mesh rendering works, is that the colour of a face is interpolated between its vertices. Now, the amount of sediment in a given column can be interpreted as sampling a location so the interpretation of colouring a vertex seems the right one. However in the depth (z) we are indicating a clear layer being **between** two depths. We need to duplicate vertices in the vertical direction in order to give ajacent layers different colours. The `explode_quad_vertices` function does exactly that:

``` {.julia #explode-vertices}
"""
    explode_quad_vertices(v)

Takes a three dimensional array representing a grid of vertices. This function duplicates these
vertices in the vertical direction, so that an amount of sediment can be given a single color.

Returns a tuple of vertices and faces (triangles), suitable for plotting with Makie's `mesh`
function.
"""
function explode_quad_vertices(v::Array{Float64,3})
    w, h, d = size(v)
    points = zeros(Float64, w, h - 1, 2, d)
    n_vertices = 2 * w * (h - 1)
    n_quads = (w - 1) * (h - 1)
    @views points[:, :, 1, :] = v[1:end, 1:end-1, :]
    @views points[:, :, 2, :] = v[1:end, 2:end, :]
    idx = reshape(1:n_vertices, w, (h - 1), 2)
    vtx1 = reshape(idx[1:end-1, :, 1], n_quads)
    vtx2 = reshape(idx[2:end, :, 1], n_quads)
    vtx3 = reshape(idx[2:end, :, 2], n_quads)
    vtx4 = reshape(idx[1:end-1, :, 2], n_quads)
    return reshape(points, n_vertices, d),
           vcat(hcat(vtx1, vtx2, vtx3), hcat(vtx1, vtx3, vtx4))
end
```

To map colours onto the returned mesh, we need to take care to duplicate each row in the input.

### Plotting Unconformities

To indicate unconformities in profile plots, we need to do two things:

1. Find those times and locations where the platform is exposed. To reduce noise we can set a minimum time span for which the platform needs to be exposed before showing up in this visualization.
2. Reduce the detected unconformities to line segments (see [the Skeleton Algorithm](../algorithms/skeleton.md)). 

``` {.julia #plot-unconformities}
"""
    plot_unconformities(ax, pack, h, minwidth; kwargs...)
    plot_unconformities(ax, pack, h, minwidth::Bool; kwargs...)
    plot_unconformities(ax, pack, h, minwidth::Int; kwargs...)

Scans the given slice for unconformities, and plots those using
Makie `linesegments`. The `minwidth` argument controls for how many time
steps the platform needs to be exposed before we plot it. For `minwidth = true`
the default width of 10 time steps is taken. The `h` argument should contain
the surface heights of the slice (see `surface_heights`).

Additional keyword arguments are forwarded to the `linesegments!` call.
"""
function plot_unconformities(ax::Axis, p::PackSlice, h, minwidth::Nothing; kwargs...)
    @info "Not plotting unconformities, got minwidth: $(minwidth)"
end

function plot_unconformities(ax::Axis, p::PackSlice, h, minwidth::Bool; kwargs...)
    if minwidth
        plot_unconformities(ax, p, h, 10; kwargs...)
    end
end

function plot_unconformities(ax::Axis, p::PackSlice, h, minwidth::Int; kwargs...)
    x = _section_positions(p)
    hiatus = skeleton(water_depth(p) .< 0.0u"m", minwidth=minwidth)

    if !isempty(hiatus[1])
        verts = [(x[pt[1]], h[pt...] |> in_units_of(u"m")) for pt in hiatus[1]]
        linesegments!(ax, vec(permutedims(verts[hiatus[2]])); kwargs...)
    end
end
```

### Plotting Co-eval Lines

``` {.julia #plot-coeval-lines}
"""
    coeval_lines!(ax::Axis, pack::PackSlice, tics::Bool)
    coeval_lines!(ax::Axis, pack::PackSlice, adm::AbstractMatrix, tics::Bool)

Plot coeval lines using default settings. If `tics` is false, nothing is done.

The `adm` argument should contain the age-depth model for the entire slice in
units of m, as computed by `surface_heights`. If the `adm` argument is left
out, it will be computed from `pack` first.

All `coeval_lines!` methods also accept a `header, data` pair instead of a
`Pack`.
"""
coeval_lines!(ax::Axis, p::PackSlice, n_tics::Bool) =
    coeval_lines!(ax, p, surface_heights(p), n_tics)

function coeval_lines!(ax::Axis, p::PackSlice, adm::AbstractMatrix{Amount}, n_tics::Bool)
    if !n_tics
        return
    else
        coeval_lines!(ax, p, adm)
    end
end

"""
    coeval_lines!(ax::Axis, pack::PackSlice, tics::Vector{Time}; kwargs...)
    coeval_lines!(ax::Axis, pack::PackSlice, adm::AbstractMatrix, tics::Vector{Time}; kwargs...)

Plot coeval lines for given times. The times in `tics` are looked up in the
time axis of `pack` to find which indices to plot.
"""
coeval_lines!(ax::Axis, p::PackSlice, tics::Vector{Time}; kwargs...) =
    coeval_lines!(ax, p, surface_heights(p), tics; kwargs...)

function coeval_lines!(ax::Axis, p::PackSlice, adm::AbstractMatrix{Amount}, tics::Vector{Time}; kwargs...)
    t = time_axis(p)
    indices::Vector{Int} = [searchsortedfirst(t, i) for i in tics]
    coeval_lines!(ax, p, adm, indices; kwargs...)
end

"""
    coeval_lines!(ax::Axis, pack::PackSlice, n_tics::Tuple{Int, Int})
    coeval_lines!(ax::Axis, pack::PackSlice, adm::AbstractMatrix, n_tics::Tuple{Int, Int})

Plot coeval lines on regular intervals given by `n_tics`. The tuple gives the number of intervals
for both minor and major tics, plotted as dotted and solid black lines respectively.
If you need more control over the aestetics of these lines, use the other `coeval_lines!` methods.
"""
coeval_lines!(ax::Axis, p::PackSlice, n_tics::Tuple{Int, Int}; kwargs...) =
    coeval_lines!(ax, p, surface_heights(p), n_tics; kwargs...)

function coeval_lines!(ax::Axis, p::PackSlice, adm::AbstractMatrix{Amount}, n_tics::Tuple{Int, Int} = (4, 8))
    n_steps = length(time_axis(p)) - 1
    n_major_tics, n_minor_tics = n_tics

    major_tics = div.(n_steps:n_steps:n_steps*n_major_tics, n_major_tics) .+ 1
    minor_tics = filter(
        t->!(t in major_tics),
        div.(n_steps:n_steps:n_steps*n_minor_tics, n_minor_tics) .+ 1) |> collect

    coeval_lines!(ax, p, adm, minor_tics, color=:black, linewidth=1, linestyle=:dot)
    coeval_lines!(ax, p, adm, major_tics, color=:black, linewidth=1, linestyle=:solid)
end

"""
    coeval_lines!(ax::Axis, pack::PackSlice, tics::Vector{Int}; kwargs...)
    coeval_lines!(ax::Axis, pack::PackSlice, adm::AbstractMatrix, tics::Vector{Int}; kwargs...)

Plot coeval lines for the given indices. The `kwargs...` are forwarded to a call to Makie's
`lines!` procedure.
"""
coeval_lines!(ax::Axis, p::PackSlice, tics::Vector{Int}; kwargs...) =
    coeval_lines!(ax, p, surface_heights(p), tics; kwargs...)

function coeval_lines!(ax::Axis, p::PackSlice, adm::AbstractMatrix{Amount}, tics::Vector{Int}; kwargs...)
    x = _section_positions(p)
    h = adm |> in_units_of(u"m")
    for t in tics
        lines!(ax, x, h[:, t]; kwargs...)
    end
end

coeval_lines!(ax::Axis, header::Header, data::DataSlice, args...; kwargs...) =
    coeval_lines!(ax, Pack(header, data), args...; kwargs...)

function plot_sealevel!(ax::Axis, p::Pack)
    sl_end = last(sea_level(header(p))) |> in_units_of(u"m")
    hlines!(ax, sl_end, color=:lightblue, linewidth=5, label="end sea level")
end
```

## Module

```{.julia file=ext/SedimentProfile.jl}
module SedimentProfile

import CarboKitten.Visualization: sediment_profile, sediment_profile!, profile_plot!, coeval_lines!, sediment_proportion!, sediment_proportion

using CarboKitten.Visualization
using CarboKitten.Utility: in_units_of
using CarboKitten.Algorithms: skeleton
using CarboKitten.Output.Storage: Header, DataSlice, Pack, PackSlice, header,
    stratigraphic_column, water_depth, surface_heights, n_facies, time_axis,
    section_axis, sea_level, bedrock

using Makie
using GeometryBasics
using Unitful

const Amount = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")


<<section-positions>>
<<explode-vertices>>
<<plot-unconformities>>
<<plot-coeval-lines>>

"""
    profile_plot!(ax, pack; color, mesh_args...)
    profile_plot!(ax, header, data_slice; color, mesh_args...)

Generic profile plot. This sets up a mesh for plotting with the Makie `mesh!`
function, plots the initial topography and the mesh by passing `mesh_args...`.

The `color` array should have the size `(n_x, n_t)`, the same as a single
facies of `production(pack)`.
"""
function profile_plot!(ax::Axis, p::PackSlice; color::AbstractArray, mesh_args...)
    x = _section_positions(p)
    h = surface_heights(p) |> in_units_of(u"m")

    verts = zeros(Float64, size(h)..., 2)
    @views verts[:, :, 1] .= x
    @views verts[:, :, 2] .= h
    v, f = explode_quad_vertices(verts)

    base = bedrock(p) |> in_units_of(u"m")
    lower_limit = minimum(base) - 20
    band!(ax, x, lower_limit, base; color=:gray, label="initial topography")
    lines!(ax, x, base; color=:black, label="initial topography")
    ylims!(ax, lower_limit + 10, nothing)
    xlims!(ax, x[1], x[end])
    ax.xlabel = "position [km]"
    ax.ylabel = "depth [m]"

    c = vec(color)
    mesh!(ax, v, f; color=vcat(c, c), mesh_args...)
end

"""
    profile_plot!(f, ax, pack; mesh_args...)
    profile_plot!(f, ax, header, data_slice; mesh_args...)

Instead of explicitely passing the color data as an array, this generates the
colors from a function `f` over the preserved sediment (stratigraphic column).
So `f` should have signature `AbstractVector -> Float` or possibly
`AbstractVector -> RGB` (whatever Makie accepts). Here the vector input has
size of the number of facies.
"""
function profile_plot!(f::F, ax::Axis, p::PackSlice; mesh_args...) where {F}
    color = f.(eachslice(stratigraphic_column(p), dims=(2, 3)))
    profile_plot!(ax, p; color=color, mesh_args...)
end

profile_plot!(ax::Axis, header::Header, data::DataSlice; kwargs...) =
    profile_plot!(ax, Pack(header, data); kwargs...)
profile_plot!(f::F, ax::Axis, header::Header, data::DataSlice; kwargs...) where {F} =
    profile_plot!(f, ax, Pack(header, data); kwargs...)

"""
    sediment_profile!(ax, pack; show_unconformities, show_coeval_lines, show_sealevel)
    sediment_profile!(ax, header, data_slice; ...)

Plot the sediment profile, choosing colour by dominant facies type (argmax).

Unconformities are shown when the sediment is subaerially exposed.
"""
function sediment_profile!(ax::Axis, p::PackSlice;
                           show_unconformities::Union{Nothing,Bool,Int} = true,
                           show_coeval_lines::Union{Bool,Tuple{Int, Int},Vector{Int},Vector{Time}} = true,
                           show_sealevel::Bool = true)
    n_f = n_facies(p)
    h = surface_heights(p)

    if show_sealevel
        plot_sealevel!(ax, p)
    end

    plot = profile_plot!(argmax, ax, p; alpha=1.0,
        colormap=cgrad(Makie.wong_colors()[1:n_f], n_f, categorical=true))

    coeval_lines!(ax, p, h, show_coeval_lines)

    minwidth = show_unconformities
    plot_unconformities(ax, p, h, minwidth; label = "unconformities",
                        color=:white, linestyle=:dash, linewidth=1)

    ax.title = "sediment profile"
    return plot
end

sediment_profile!(ax::Axis, header::Header, data::DataSlice; kwargs...) =
    sediment_profile!(ax, Pack(header, data); kwargs...)

"""
    sediment_profile(pack; show_unconformities=true)
    sediment_profile(header, data_slice; show_unconformities=true)

Plot the sediment profile from a slice. Dominant facies colour is chosen
by `argmax`. By default unconformities are shown using dashed white lines.
"""
function sediment_profile(p::PackSlice;
                          show_unconformities::Union{Bool,Int,Nothing} = true)
    fig = Figure(size=(1000, 600))
    ax = Axis(fig[1, 1])
    sediment_profile!(ax, p; show_unconformities=show_unconformities)
    return fig
end

sediment_profile(header::Header, data::DataSlice; kwargs...) =
    sediment_profile(Pack(header, data); kwargs...)

"""
    sediment_proportion!(ax, pack, facies_index; colorrange=(0,1), colormap=:viridis)
    sediment_proportion!(ax, header, data, facies_index; ...)

Plot the proportion of `facies_index` relative to total sediment at each
location and time step, using the same stratigraphic mesh as `sediment_profile!`.

Returns the `mesh!` plot object (for attaching a `Colorbar`).
"""
function sediment_proportion!(ax::Axis, p::PackSlice, facies_index::Int;
                               colorrange::Tuple = (0.0, 1.0),
                               colormap = :viridis)
    n_f = n_facies(p)
    @assert 1 <= facies_index <= n_f "facies_index $(facies_index) out of range 1:$(n_f)"

    proportion = map(eachslice(stratigraphic_column(p), dims=(2, 3))) do col
        total = sum(col)
        total > zero(total) ? col[facies_index] / total : 0.0
    end

    plot = profile_plot!(ax, p; color=proportion,
        colorrange=colorrange, colormap=colormap)
    ax.title = "proportion of facies $(facies_index)"
    return plot
end

sediment_proportion!(ax::Axis, header::Header, data::DataSlice, facies_index::Int; kwargs...) =
    sediment_proportion!(ax, Pack(header, data), facies_index; kwargs...)

"""
    sediment_proportion(pack, facies_index; kwargs...)
    sediment_proportion(header, data, facies_index; kwargs...)

Standalone proportion figure. Shows the fraction of `facies_index` relative to
total sediment using the same mesh geometry as `sediment_profile`. A `Colorbar`
is added automatically.

See `sediment_proportion!` for keyword arguments.
"""
function sediment_proportion(p::PackSlice, facies_index::Int;
                              colorrange::Tuple = (0.0, 1.0),
                              colormap = :viridis)
    fig = Figure(size=(1000, 600))
    ax  = Axis(fig[1, 1])
    plot = sediment_proportion!(ax, p, facies_index;
        colorrange=colorrange, colormap=colormap)
    Colorbar(fig[1, 2], plot; label="facies $(facies_index) proportion")
    return fig
end

sediment_proportion(header::Header, data::DataSlice, facies_index::Int; kwargs...) =
    sediment_proportion(Pack(header, data), facies_index; kwargs...)

end  # module
```
