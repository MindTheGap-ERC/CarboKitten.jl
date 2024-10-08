# ~/~ begin <<docs/src/ca-with-production.md#ext/VisualizationExt.jl>>[init]
module VisualizationExt

import CarboKitten.Visualization: plot_facies_production, plot_crosssection, sediment_profile

using CarboKitten
using CarboKitten.Visualization
using CarboKitten.Burgess2013: production_rate
using CarboKitten.Utility: in_units_of

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

struct Axes
	x::Vector{Length}
	y::Vector{Length}
	t::Vector{Time}
end

struct Header
	axes::Axes
	Δt::Time
	time_steps::Int
	bedrock_elevation::Matrix{Amount}
	sea_level::Vector{Length}
	subsidence_rate::Rate
end

struct Data
	disintegration::Array{Amount, 4}
	production::Array{Amount, 4}
	deposition::Array{Amount, 4}
	sediment_elevation::Array{Amount, 3}
end

struct DataSlice
	disintegration::Array{Amount, 3}
	production::Array{Amount, 3}
	deposition::Array{Amount, 3}
	sediment_elevation::Array{Amount, 2}
end

function read_header(fid)
	attrs = HDF5.attributes(fid["input"])

	axes = Axes(
		fid["input/x"][]*u"m",
		fid["input/y"][]*u"m",
		fid["input/t"][]*u"Myr")

	return Header(
		axes,
		attrs["delta_t"][]*u"Myr",
		attrs["time_steps"][],
		fid["input/bedrock_elevation"][]*u"m",
		fid["input/sea_level"][]*u"m",
		attrs["subsidence_rate"][]*u"m/Myr")
end

function read_data(filename)
	h5open(filename) do fid
		header = read_header(fid)
		data = Data(
			fid["disintegration"][]*u"m",
			fid["production"][]*u"m",
			fid["deposition"][]*u"m",
			fid["sediment_height"][]*u"m")
		header, data
	end
end

function read_slice(filename, slice...)
	h5open(filename) do fid
		header = read_header(fid)
		data = DataSlice(
			fid["disintegration"][slice...]*u"m",
			fid["production"][slice...]*u"m",
			fid["deposition"][slice...]*u"m",
			fid["sediment_height"][slice[2:end]...]*u"m")
		header, data
	end
end

elevation(h::Header, d::Data) = let bl = h.bedrock_elevation[:,:,na],
									sr = h.axes.t[end] * h.subsidence_rate
	cat(bl, bl .+ d.sediment_elevation; dims=3) .- sr
end

elevation(h::Header, d::DataSlice, y) = let bl = h.bedrock_elevation[:,y,na],
	       								 sr = h.axes.t[end] * h.subsidence_rate
	cat(bl, bl .+ d.sediment_elevation; dims=2) .- sr
end

colormax(d::Data) = getindex.(argmax(d.deposition; dims=1)[1,:,:,:], 1)
colormax(d::DataSlice) = getindex.(argmax(d.deposition; dims=1)[1,:,:], 1)

"""
    explode_quad_vertices(v)

Takes a three dimensional array representing a grid of vertices. This function duplicates these
vertices in the vertical direction, so that an amount of sediment can be given a single color.

Returns a tuple of vertices and faces (triangles), suitable for plotting with Makie's `mesh`
function.
"""
function explode_quad_vertices(v::Array{Float64, 3})
	w, h, d = size(v)
	points = zeros(Float64, w, h-1, 2, d)
	n_vertices = 2 * w * (h-1)
	n_quads = (w - 1) * (h - 1)
	@views points[:, :, 1, :] = v[1:end, 1:end-1, :]
	@views points[:, :, 2, :] = v[1:end, 2:end, :]
	idx = reshape(1:n_vertices, w, (h-1), 2)
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
	dxs = CartesianIndices(ntuple(_->3, dim)) .|> (x -> x - CartesianIndex(ntuple(_->2, dim))) |> filter(x->x!=CartesianIndex(ntuple(_->0, dim)...))
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
	return out, group-1
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

	water_depth = ξ .- (header.subsidence_rate .* (header.axes.t .- header.axes.t[end]) .+ header.sea_level)[na, :]
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
	fig = Figure(size=(800,600))
	ax = Axis(fig[1,1])
	sediment_profile!(ax, filename, y)
	return fig
end

function plot_facies_production(input; loc = nothing)
 fig, loc = isnothing(loc) ? let fig = Figure(); (fig, fig[1, 1]) end : (nothing, loc)
 ax = Axis(loc, title="production at $(sprint(show, input.insolation; context=:fancy_exponent=>true))", xlabel="production (m/Myr)", ylabel="depth (m)", yreversed=true)
 for f in input.facies
  depth = (0.1:0.1:50.0)u"m"
  prod = [production_rate(input.insolation, f, d) for d in depth]
  lines!(ax, prod / u"m/Myr", depth / u"m")

 end
 fig
end

function plot_crosssection(pos, datafile)
    # x: 1-d array with x-coordinates
    # t: 1-d array with time-coordinates (n_steps + 1)
    # h[x, t]: height fn, monotonic increasing in time
    # p[x, facies, t]: production rate
    # taken at y = y_max / 2, h[x, 1] is initial height
    n_facies, x, t, h, p = h5open(datafile, "r") do fid
        attr = HDF5.attributes(fid["input"])
        Δt = attr["delta_t"][]
        subsidence_rate = attr["subsidence_rate"][]
        t_end = fid["input/t"][end-1]
        total_subsidence = subsidence_rate * t_end
        total_sediment = sum(fid["sediment"][]; dims=3)
        initial_height = fid["input/height"][]
        center = div(size(total_sediment)[1], 2)
        elevation = cumsum(total_sediment; dims=4)[:, center, 1, :] .* Δt .- initial_height .- total_subsidence
        t = fid["input/t"][]
        n_facies = size(fid["sediment"])[3]

        return n_facies,
        fid["input/x"][],
        [t; Δt * attr["time_steps"][]],
        hcat(.-initial_height .- total_subsidence, elevation),
        fid["sediment"][:, center, :, :]
    end

    pts = vec(Point{2,Float64}.(x, h[:, 2:end]))
    c = vec(argmax(p; dims=2)[:, 1, :] .|> (c -> c[2]))
    rect = Rect2(0.0, 0.0, 1.0, 1.0)
    m_tmp = GeometryBasics.mesh(Tesselation(rect, (100, 1000)))
    m = GeometryBasics.Mesh(pts, faces(m_tmp))

    # pts = vec(Point{2,Float64}.(x, h))
    # c = argmax(p; dims=2)[:,1,:] .|> (c -> c[2])
    # w = size(x)[1]

    # face(idx) = let k = idx[1] + idx[2]*w
    #     TriangleFace(k, k+1, k+1+w), TriangleFace(k+1+w, k+w, k)
    # end

    ax = Axis(pos, xlabel="location", ylabel="depth", limits=((-12, x[end]), nothing))
    # for f in 1:n_facies
    #     locs = CartesianIndices((size(x)[1], size(t)[1] - 1))[c .== f]
    #     triangles = collect(Iterators.flatten(face.(locs)))
    #     m = GeometryBasics.Mesh(pts, triangles)
    #     mesh!(ax, m)
    # end

    mesh!(ax, m, color=c, alpha=0.7)
    for idx in [1, 501, 1001]
        lines!(ax, x, h[:, idx], color=:black)
        text!(ax, -2.0, h[1, idx]; text="$(t[idx]) Myr", align=(:right, :center))
    end
    for idx in [250, 750]
        lines!(ax, x, h[:, idx], color=:black, linewidth=0.5)
    end
end

end
# ~/~ end
