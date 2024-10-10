### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 51490338-86e9-11ef-106c-5bede5bc2b64
using Pkg; Pkg.activate("../workenv")

# ╔═╡ b080be45-7979-4214-b64e-a0fbed0ec93b
using Revise

# ╔═╡ a46f4043-f471-41d5-9c81-2396031f9a89
using GLMakie

# ╔═╡ bc4fd4d7-5b16-4efa-9e9f-a1ec96ebc9af
using Unitful

# ╔═╡ 4bdc4267-5c95-4d9a-ae61-2a98ddd6717b
using CarboKitten.Export: read_slice, Header, DataSlice

# ╔═╡ a67445cf-1947-4dd1-a75a-e8d1b00eaccf
using CarboKitten.Visualization

# ╔═╡ 359483a3-f2f9-4470-ba11-730fc9bf8b46
md"""
In the case of the BS92 model, the idea of marking sedimentation gaps in the sediment profile using connected components doesn't work. We get the following sediment curves.
"""

# ╔═╡ 5db8cfd7-7d45-4d77-abf8-ac0e1270d04a
let
	header, data = read_slice("../data/output/bs92.h5", :, 1)
	fig = Figure()
	ax = Axis(fig[1,1])
	for (i, t) in enumerate(header.axes.t)
		η = header.bedrock_elevation .+ data.sediment_elevation[:, i]
		lines!(ax, header.axes.x, η[:,1]; color=Makie.wong_colors()[1])
	end
	fig
end

# ╔═╡ 49404961-56c5-455e-810e-afc59ede0148
md"""Meanwhile, the `sediment_profile` routine creates something non-sensical"""

# ╔═╡ dbce83bb-066d-41ee-a60d-5455c33b337f
let
	header, data = read_slice("../data/output/bs92.h5", :, 1)
	sediment_profile(header, data)
end

# ╔═╡ a7f7077f-751d-4377-b43b-c9520fa4f768
md"""In the Wheeler diagram we see that all the gaps are connected. So we need something more advanced to render the gaps."""

# ╔═╡ 2ea5a1c8-2bf5-4f59-8bb7-0f2d06370812
let
	header, data = read_slice("../data/output/bs92.h5", :, 1)
	wheeler_diagram(header, data; smooth_size=(1,1), range=(-5000.0u"m/Myr", 5000.0u"m/Myr"))
end

# ╔═╡ eddc0e5e-65f0-4420-905a-59660d4aec4f
const na = [CartesianIndex()]

# ╔═╡ e91e092e-bb97-4d93-8702-a6bb650a094f
elevation(h::Header, d::DataSlice) =
    let bl = h.bedrock_elevation[d.slice..., na],
        sr = h.axes.t[end] * h.subsidence_rate

        bl .+ d.sediment_elevation .- sr
    end

# ╔═╡ c46c6eea-eaf2-4d79-b478-e7861ef3f07f
water_depth(h::Header, d::DataSlice) = elevation(h, d) .- (h.subsidence_rate.*(h.axes.t.-h.axes.t[end]).+h.sea_level)[na, :]

# ╔═╡ 8fe2d845-14c0-4154-982d-412dd9e0140a
let
	header, data = read_slice("../data/output/bs92.h5", :, 1)
	w = water_depth(header, data)

	fig = Figure()
	ax = Axis(fig[2,1])
	
	hm = heatmap!(ax, header.axes.x/u"km" .|> NoUnits, header.axes.t/u"Myr" .|> NoUnits, w/u"m" .|> NoUnits)
	contour!(ax, header.axes.x/u"km" .|> NoUnits, header.axes.t/u"Myr" .|> NoUnits, w/u"m" .|> NoUnits; levels=[0.0], color=:black)
	Colorbar(fig[1, 1], hm; vertical=false, label="water depth [m]")
	fig
end

# ╔═╡ 1b4c0aab-4a76-4255-8fca-173202a9c068
wd = water_depth(read_slice("../data/output/bs92.h5", :, 1)...)

# ╔═╡ 35c58077-5237-429f-99a2-9d974297ad22
dry = wd .> 0.0u"m"

# ╔═╡ 8d98914b-c8f6-4937-8d62-284fbf82c4de
heatmap(dry)

# ╔═╡ 579e7112-4cc2-4488-a993-5cae889cf84f
struct RangeFinder
	v::AbstractVector{Bool}
end

# ╔═╡ 9d654482-3eb2-41c7-b57d-f47abb1605af
Base.iterate(r::RangeFinder) = iterate(r, 1)

# ╔═╡ 26021a8e-8e0c-4fba-a3f0-247e1fbc99d9
Base.eltype(r::RangeFinder) = UnitRange{Int}

# ╔═╡ 657c153a-b8e7-478f-9120-a14a95fae771
Base.IteratorSize(::Type{RangeFinder}) = Base.SizeUnknown()

# ╔═╡ c0dc5de7-70f7-4384-b730-9aa19aef1301
RangeFinder(dry[50,:]) |> collect

# ╔═╡ 1c0867c2-ea0c-46e8-8ab3-24e98d764ee0
struct TagVectors{T}
	vectors::T
end

# ╔═╡ 29717b91-e16a-460b-8255-d15b2cf7aefc
Base.size(r::TagVectors{T}) where T = size(r.vectors)

# ╔═╡ 32f63a12-b5ad-413d-86a3-b71642ce1971
Base.length(r::TagVectors{T}) where T = length(r.vectors)

# ╔═╡ 031fcdfb-7191-4eb6-af7a-f9c5c6864b25
Base.eltype(::Type{TagVectors{T}}) where T = Any
# This should be:
# Iterators.Zip{Tuple{UnitRange{Int},Vector{eltype(eltype(T))}}}
# But that doesn't work

# ╔═╡ 2fd0599c-b87c-47b5-8bd6-bb51eeb0af7e
begin
	struct Pairs{T}
		v::T
	end

	function Base.iterate(p::Pairs{T}) where T
		x = iterate(p.v)
		isnothing(x) && return nothing
		return iterate(p, x)
	end

	function Base.iterate(p::Pairs{T}, st) where T
		(prev, it) = st
		x = iterate(p.v, it)
		isnothing(x) && return nothing
		curr, nit = x
		return (prev, curr), (curr, nit)
	end

	Base.eltype(::Type{Pairs{T}}) where T = NTuple{2, eltype(T)}
	Base.IteratorSize(::Type{Pairs{T}}) where T = Base.IteratorSize(T)
	Base.size(p::Pairs{T}) where T = length(p)
	Base.length(p::Pairs{T}) where T = length(p.v) - 1
end

# ╔═╡ 1f02d0ac-6703-4951-86d1-1894551edca8
function Base.iterate(r::RangeFinder, i::Union{Int, Nothing})
	isnothing(i) && return nothing
	a = findnext(r.v, i)
	isnothing(a) && return nothing
	b = findnext(!, r.v, a)
	isnothing(b) && return (a:length(r.v)), nothing
	return (a:b-1), b
end

# ╔═╡ 20fe3b99-77c1-4ff5-93ac-8bd1fe55b7d6
function Base.iterate(tv::TagVectors{T}) where T
	tag = 1
	x = iterate(tv.vectors)
	isnothing(x) && return nothing
	(v_, nit) = x
	v = collect(v_)
	n = length(v)
	return zip(tag:tag+n-1, v) |> collect, (tag+n, nit)
end

# ╔═╡ 6db7d714-1982-419e-98b9-ff86b2f0a127
function Base.iterate(tv::TagVectors{T}, st::Tuple{Int,U}) where {T, U}
	(tag, it) = st
	isnothing(it) && return nothing
	x = iterate(tv.vectors, it)
	isnothing(x) && return nothing
	(v_, nit) = x
	v = collect(v_)
	n = length(v)
	return zip(tag:tag+n-1, v) |> collect, (tag+n, nit)
end

# ╔═╡ 96edc0b3-1b96-430d-9821-a00ffbe34ad1
Base.IteratorSize(::Type{TagVectors{T}}) where T = Base.IteratorSize(T)

# ╔═╡ a9136d78-35be-4d46-999f-8e9238abdf79
Pairs([1, 2, 3]) |> collect

# ╔═╡ 1ab6b19b-6684-4aae-80b5-f5c8107ae25f
pairs(it) = zip(it, Iterators.drop(it, 1))

# ╔═╡ 575243e0-6f89-4f1f-b278-d39bf1715b9b
pairs([1, 2, 3]) |> collect

# ╔═╡ 1f69a952-3899-410d-80a1-2955b3d1d99d
const Vertex = Tuple{Int, UnitRange{Int}}

# ╔═╡ 1116dad5-9c32-4d81-b8eb-2c5eb6fe06fe
isempty((1:10) ∩ (13:15))

# ╔═╡ b32736e5-bc0b-463a-8eb5-5099df519217
edge(a::Vertex, b::Vertex) = isempty(a[2] ∩ b[2]) ? nothing : (a[1], b[1])

# ╔═╡ d63b489e-04f9-4a02-a78a-f530acb964b3
edges(a, b) = Iterators.filter(!isnothing, Iterators.map(splat(edge), Iterators.product(a, b)))

# ╔═╡ 40be425e-7912-4d7f-ab1d-cb1ed2e066d2
middle(a::UnitRange{Int}) = (a.start + a.stop) ÷ 2

# ╔═╡ 83504127-77aa-4d83-b00b-9301faa830a8
(vtx, egs) = let
	wd = water_depth(read_slice("../data/output/bs92.h5", :, 1)...)
	dry = wd .> 0.0u"m"
	vts = Iterators.map(RangeFinder, eachrow(dry))
	egs = Iterators.flatten(Iterators.map(splat(edges), Pairs(vts |> TagVectors))) |> collect
	vtx = Iterators.flatten(((i, middle(v)) for v in vs) for (i, vs) in enumerate(vts)) |> collect
	vtx, reshape(reinterpret(Int, egs), (2,:))'
end

# ╔═╡ 3b08d18e-69af-44e5-97dc-13e35dafba28
let
	fig, ax = heatmap(dry)
	linesegments!(ax, vec(permutedims(vtx[egs])), color=:black)
	fig
end

# ╔═╡ e8756327-6533-4cf4-acfb-58019ad74c7d


# ╔═╡ Cell order:
# ╠═51490338-86e9-11ef-106c-5bede5bc2b64
# ╠═b080be45-7979-4214-b64e-a0fbed0ec93b
# ╠═a46f4043-f471-41d5-9c81-2396031f9a89
# ╠═bc4fd4d7-5b16-4efa-9e9f-a1ec96ebc9af
# ╠═4bdc4267-5c95-4d9a-ae61-2a98ddd6717b
# ╠═a67445cf-1947-4dd1-a75a-e8d1b00eaccf
# ╟─359483a3-f2f9-4470-ba11-730fc9bf8b46
# ╟─5db8cfd7-7d45-4d77-abf8-ac0e1270d04a
# ╟─49404961-56c5-455e-810e-afc59ede0148
# ╠═dbce83bb-066d-41ee-a60d-5455c33b337f
# ╟─a7f7077f-751d-4377-b43b-c9520fa4f768
# ╟─2ea5a1c8-2bf5-4f59-8bb7-0f2d06370812
# ╠═eddc0e5e-65f0-4420-905a-59660d4aec4f
# ╠═e91e092e-bb97-4d93-8702-a6bb650a094f
# ╠═c46c6eea-eaf2-4d79-b478-e7861ef3f07f
# ╠═8fe2d845-14c0-4154-982d-412dd9e0140a
# ╠═1b4c0aab-4a76-4255-8fca-173202a9c068
# ╠═35c58077-5237-429f-99a2-9d974297ad22
# ╠═8d98914b-c8f6-4937-8d62-284fbf82c4de
# ╠═579e7112-4cc2-4488-a993-5cae889cf84f
# ╠═9d654482-3eb2-41c7-b57d-f47abb1605af
# ╠═1f02d0ac-6703-4951-86d1-1894551edca8
# ╠═26021a8e-8e0c-4fba-a3f0-247e1fbc99d9
# ╠═657c153a-b8e7-478f-9120-a14a95fae771
# ╠═c0dc5de7-70f7-4384-b730-9aa19aef1301
# ╠═1c0867c2-ea0c-46e8-8ab3-24e98d764ee0
# ╠═20fe3b99-77c1-4ff5-93ac-8bd1fe55b7d6
# ╠═6db7d714-1982-419e-98b9-ff86b2f0a127
# ╠═96edc0b3-1b96-430d-9821-a00ffbe34ad1
# ╠═29717b91-e16a-460b-8255-d15b2cf7aefc
# ╠═32f63a12-b5ad-413d-86a3-b71642ce1971
# ╠═031fcdfb-7191-4eb6-af7a-f9c5c6864b25
# ╠═2fd0599c-b87c-47b5-8bd6-bb51eeb0af7e
# ╠═a9136d78-35be-4d46-999f-8e9238abdf79
# ╠═1ab6b19b-6684-4aae-80b5-f5c8107ae25f
# ╠═575243e0-6f89-4f1f-b278-d39bf1715b9b
# ╠═1f69a952-3899-410d-80a1-2955b3d1d99d
# ╠═1116dad5-9c32-4d81-b8eb-2c5eb6fe06fe
# ╠═b32736e5-bc0b-463a-8eb5-5099df519217
# ╠═d63b489e-04f9-4a02-a78a-f530acb964b3
# ╠═40be425e-7912-4d7f-ab1d-cb1ed2e066d2
# ╠═83504127-77aa-4d83-b00b-9301faa830a8
# ╠═3b08d18e-69af-44e5-97dc-13e35dafba28
# ╠═e8756327-6533-4cf4-acfb-58019ad74c7d
