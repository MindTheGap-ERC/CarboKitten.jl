# Visualization

The visualization of CarboKitten output is implemented in a Julia [package extension](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)). This is done so that `CarboKitten.jl` itself doesn't have to depend on `Makie.jl` (our main visualization tool), which has a large transient dependency stack. To make the `Visualization` extension of CarboKitten available, make sure to activate a Julia project where `Makie` is installed.

## Makie primer

`Makie.jl` is a visualization package that creates exceptionally good looking (publication quality) plots in both 2D and 3D. There are three back-ends for Makie:

- `CairoMakie` for publication quality vector graphics, writing to `SVG`, `PDF` or `PNG`.
- `GLMakie` has better run-time performance than `CairoMakie`, especially when dealing with larger datasets and/or 3D visualizations. However, `GLMakie` can only produce rasterized images, so `PNG`, `JPEG` or directly to screen for interactive use.
- `WGLMakie` for online publication using WebGL. If you want interactive plots, like 3D plots that you can rotate in the browser, this is the one to use. Fair warning: this is also the least stable back-end for Makie.

To work with Makie, you need to import one of the three back-end packages. In general, every plot available in Makie has two variants. One is a direct function for plotting:

```julia
using CairoMakie

x = randn(10)
y = randn(10)

scatter(x, y)
```

The other requires a bit more prep, but gives you more control.

```julia
fig = Figure()
ax = Axis(fig[1,1])
scatter!(ax, x, y)
```

Here, we create a figure explicitly, then create a new set of axes somewhere on the grid in the figure, and then plot on that set of axes. The plotting functions accepting an `Axis` argument actually modify an existing context, which is why these functions always end with an exclamation mark, in this case `scatter!`.

If you like to know more about Makie, their ["Getting started"](https://docs.makie.org/stable/tutorials/getting-started) is a good place to start.

## Project Extension

The Project Extension requires a front-end where the available methods are exposed.

``` {.julia file=src/Visualization.jl}
module Visualization
export sediment_profile!, sediment_profile, wheeler_diagram!, wheeler_diagram, production_curve!, production_curve

print_instructions() = print("This is an extension and only becomes available when you import {Cairo,GL,WGL}Makie before using this.\n")

sediment_profile!(args...) = print_instructions()
sediment_profile(args...) = print_instructions()
wheeler_diagram!(args...) = print_instructions()
wheeler_diagram(args...) = print_instructions()
production_curve(args...) = print_instructions()
production_curve!(args...) = print_instructions()
end  # module
```

## Wheeler diagram

``` {.julia file=ext/WheelerDiagram.jl}
module WheelerDiagram

import CarboKitten.Visualization: wheeler_diagram
using CarboKitten.Export: Header, Data, DataSlice, read_data, read_slice
using CarboKitten.Utility: in_units_of
using Makie
using Unitful
using CarboKitten.BoundaryTrait
using CarboKitten.Stencil: convolution


const na = [CartesianIndex()]

elevation(h::Header, d::DataSlice) =
    let bl = h.bedrock_elevation[d.slice..., na],
        sr = h.axes.t[end] * h.subsidence_rate

        cat(bl, bl .+ d.sediment_elevation; dims=2) .- sr
    end

water_depth(header::Header, data::DataSlice) =
    let h = elevation(header, data),
        s = header.subsidence_rate .* (header.axes.t .- header.axes.t[end]),
        l = header.sea_level

        h .- (s .+ l)[na, :]
    end

const Rate = typeof(1.0u"m/Myr")

function sediment_accumulation!(ax::Axis, header::Header, data::DataSlice;
                                smooth_size::NTuple{2,Int} = (3, 11),
                                colormap = Reverse(:curl),
                                range::NTuple{2, Rate} = (-100.0u"m/Myr", 100.0u"m/Myr"))
	magnitude = sum(data.deposition .- data.disintegration; dims=1)[1, :, :] ./ header.Δt
    blur = convolution(Shelf, ones(Float64, smooth_size...) ./ *(smooth_size...))
	wd = zeros(Float64, length(header.axes.x), length(header.axes.t))
	blur(water_depth(header, data) / u"m", wd)
	mag = zeros(Float64, length(header.axes.x), length(header.axes.t) - 1)
	blur(magnitude / u"m/Myr", mag)

	ax.ylabel = "time [Myr]"
	ax.xlabel = "position [km]"

	sa = heatmap!(ax, header.axes.x / u"km", header.axes.t / u"Myr", mag;
        colormap=colormap, colorrange=range./u"m/Myr")
	contour!(ax, header.axes.x / u"km", header.axes.t / u"Myr", wd;
        levels=[0], color=:red, linewidth=2, linestyle=:dash)
    return sa
end

function dominant_facies!(ax::Axis, header::Header, data::DataSlice;
                          smooth_size::NTuple{2,Int} = (3, 11),
                          colors = Makie.wong_colors())
    n_facies = size(data.production)[1]
    colormax(d) = getindex.(argmax(d; dims=1)[1, :, :], 1)

	dominant_facies = colormax(data.deposition)
	blur = convolution(Shelf, ones(Float64, 3, 11) ./ 33.0)
	wd = zeros(Float64, length(header.axes.x), length(header.axes.t))
	blur(water_depth(header, data) / u"m", wd)

	ax.ylabel = "time [Myr]"
	ax.xlabel = "position [km]"

	ft = heatmap!(ax, header.axes.x / u"km", header.axes.t / u"Myr", dominant_facies;
	    colormap = cgrad(colors[1:n_facies], n_facies, categorical=true),
        colorrange=(0.5, n_facies + 0.5))
	contourf!(ax, header.axes.x / u"km", header.axes.t / u"Myr", wd;
        levels=[0.0, 100.0], colormap=Reverse(:grays))
	contour!(ax, header.axes.x / u"km", header.axes.t / u"Myr", wd;
        levels=[0], color=:black, linewidth=2)
    return ft	
end

function wheeler_diagram(header::Header, data::DataSlice)
	fig = Figure(size=(1000, 600))
	ax1 = Axis(fig[2,1])
	ax2 = Axis(fig[2,2])

	linkyaxes!(ax1, ax2)

    sa = sediment_accumulation!(ax1, header, data)
    ft = dominant_facies!(ax2, header, data)
    ax2.ylabel = ""

	Colorbar(fig[1,1], sa; vertical=false, label="sediment accumulation [m/Myr]")
	Colorbar(fig[1,2], ft; vertical=false, ticks=1:3, label="dominant facies")
	fig
end

end
```

## Production curve

``` {.julia file=ext/ProductionCurve.jl}
module ProductionCurve

using Makie
using Unitful

import CarboKitten.Visualization: production_curve!, production_curve
using CarboKitten.Burgess2013: production_rate

function production_curve!(ax, input; loc=nothing)
end

function production_curve(input; loc=nothing)
    fig, loc = isnothing(loc) ? let fig = Figure()
        (fig, fig[1, 1])
    end : (nothing, loc)
    ax = Axis(loc, title="production at $(sprint(show, input.insolation; context=:fancy_exponent=>true))",
              xlabel="production (m/Myr)", ylabel="depth (m)", yreversed=true)
    for f in input.facies
        depth = (0.1:0.1:50.0)u"m"
        prod = [production_rate(input.insolation, f, d) for d in depth]
        lines!(ax, prod / u"m/Myr", depth / u"m")

    end
    fig
end

end
```

## Sediment profile

``` {.julia file=ext/VisualizationExt.jl}
module VisualizationExt

include("WheelerDiagram.jl")
include("ProductionCurve.jl")

import CarboKitten.Visualization: sediment_profile, sediment_profile!

using CarboKitten
using CarboKitten.Visualization
using CarboKitten.Burgess2013: production_rate
using CarboKitten.Utility: in_units_of
using CarboKitten.Export: Header, Data, DataSlice, read_data, read_slice

using HDF5
using Makie
using GeometryBasics
using Unitful

using Statistics: mean

const Rate = typeof(1.0u"m/Myr")
const Amount = typeof(1.0u"m")
const Length = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")

const na = [CartesianIndex()]

elevation(h::Header, d::Data) =
    let bl = h.bedrock_elevation[:, :, na],
        sr = h.axes.t[end] * h.subsidence_rate

        cat(bl, bl .+ d.sediment_elevation; dims=3) .- sr
    end

elevation(h::Header, d::DataSlice) =
    let bl = h.bedrock_elevation[d.slice..., na],
        sr = h.axes.t[end] * h.subsidence_rate

        cat(bl, bl .+ d.sediment_elevation; dims=2) .- sr
    end

colormax(d::Data) = getindex.(argmax(d.deposition; dims=1)[1, :, :, :], 1)
colormax(d::DataSlice) = getindex.(argmax(d.deposition; dims=1)[1, :, :], 1)

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

"""
    bean_counter(mask)

Given a mask (array of binary values), performs a floodfill on all connected components,
giving each an integer identifier.

Returns the array of integers identifying each group and the number of groups.
"""
function bean_counter(mask::BitArray{dim}) where {dim}
    visited = BitArray{dim}(undef, size(mask)...)
    visited .= false
    out = zeros(Int, size(mask)...)
    dxs = CartesianIndices(ntuple(_ -> 3, dim)) .|> (x -> x - CartesianIndex(ntuple(_ -> 2, dim))) |> filter(x -> x != CartesianIndex(ntuple(_ -> 0, dim)...))
    group = 1

    for idx in CartesianIndices(mask)
        visited[idx] && continue
        visited[idx] = true
        mask[idx] || continue
        out[idx] = group

        stack = idx .+ dxs
        while !isempty(stack)
            jdx = pop!(stack)
            checkbounds(Bool, mask, jdx) || continue
            visited[jdx] && continue
            visited[jdx] = true
            mask[jdx] || continue
            out[jdx] = group
            append!(stack, jdx .+ dxs)
        end
        group += 1
    end
    return out, group - 1
end

function sediment_profile!(ax::Axis, filename::AbstractString, y::Int)
    header, data = read_slice(filename, :, y)
    sediment_profile!(ax, header, data)
end

function sediment_profile!(ax::Axis, header::Header, data::DataSlice)
    x = header.axes.x |> in_units_of(u"km")
    t = header.axes.t |> in_units_of(u"Myr")
    n_facies = size(data.production)[1]
    ξ = elevation(header, data)  # |> in_units_of(u"m")

    verts = zeros(Float64, length(x), length(t), 2)
    @views verts[:, :, 1] .= x
    @views verts[:, :, 2] .= ξ |> in_units_of(u"m")
    v, f = explode_quad_vertices(verts)

    water_depth = ξ .- (header.subsidence_rate.*(header.axes.t.-header.axes.t[end]).+header.sea_level)[na, :]
    gaps, n_gaps = bean_counter(water_depth .> 0u"m")

    total_subsidence = header.subsidence_rate * header.axes.t[end]
    bedrock = (header.bedrock_elevation[data.slice...] .- total_subsidence) |> in_units_of(u"m")
    lower_limit = minimum(bedrock) - 20
    band!(ax, x, lower_limit, bedrock; color=:gray)
    lines!(ax, x, bedrock; color=:black)
    ylims!(ax, lower_limit + 10, nothing)
    xlims!(ax, x[1], x[end])
    ax.xlabel = "position [km]"
    ax.ylabel = "height [m]"
    ax.title = "sediment profile"

    c = reshape(colormax(data)[:, :], length(x) * (length(t) - 1))
    mesh!(ax, v, f, color=vcat(c, c), alpha=1.0, colormap=cgrad(Makie.wong_colors()[1:n_facies], n_facies, categorical=true))

    for g = 1:n_gaps
        size = sum(gaps .== g)
        if size < 1000
            continue
        end
        # compute the mean z-value for a gap
        gap = mean.(skipmissing.(eachslice(CartesianIndices(ξ) .|> (i -> gaps[i] == g ? ξ[i] : missing), dims=(1,))))
        lines!(ax, x, gap |> in_units_of(u"m"), color=:white, linewidth=2, linestyle=:dash)
    end
end

function sediment_profile(header::Header, data_slice::DataSlice)
    fig = Figure(size=(1000, 600))
    ax = Axis(fig[1, 1])
    sediment_profile!(ax, header, data_slice)
    return fig
end

function sediment_profile(filename, y)
    fig = Figure(size=(1000, 600))
    ax = Axis(fig[1, 1])
    sediment_profile!(ax, filename, y)
    return fig
end

end
```
