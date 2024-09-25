# ~/~ begin <<docs/src/visualization.md#ext/VisualizationExt.jl>>[init]
module VisualizationExt

import CarboKitten.Visualization: plot_facies_production, plot_crosssection, sediment_profile, sediment_profile!

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

elevation(h::Header, d::DataSlice, y) =
    let bl = h.bedrock_elevation[:, y, na],
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

function sediment_profile!(ax::Axis, filename, y)
    header, data = read_slice(filename, :, :, y, :)
    x = header.axes.x |> in_units_of(u"km")
    t = header.axes.t |> in_units_of(u"Myr")
    ξ = elevation(header, data, y)  # |> in_units_of(u"m")

    verts = zeros(Float64, length(x), length(t), 2)
    @views verts[:, :, 1] .= x
    @views verts[:, :, 2] .= ξ |> in_units_of(u"m")
    v, f = explode_quad_vertices(verts)

    water_depth = ξ .- (header.subsidence_rate.*(header.axes.t.-header.axes.t[end]).+header.sea_level)[na, :]
    gaps, n_gaps = bean_counter(water_depth .> 0u"m")

    c = reshape(colormax(data)[:, :], length(x) * (length(t) - 1))
    mesh!(ax, v, f, color=vcat(c, c), alpha=1.0)

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

function sediment_profile(filename, y)
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1])
    sediment_profile!(ax, filename, y)
    return fig
end

function plot_facies_production(input; loc=nothing)
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
# ~/~ end
