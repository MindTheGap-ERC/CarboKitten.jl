### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 4aaa0aee-393d-11ef-1357-739a4dfddae9
using Pkg ; Pkg.activate("../workenv")

# ╔═╡ 30b991e7-5034-46c1-a46b-990f8055e4c8
using HDF5

# ╔═╡ 3e9d6678-1625-46c5-91c2-85f231d668a5
using Unitful

# ╔═╡ 3d13974e-c548-4303-bdba-c09f241fbe59
using GLMakie

# ╔═╡ c45863cd-1f26-4039-bd1b-37b4e7c89ec6
using CarboKitten.Model.ALCAPS

# ╔═╡ 59fcb379-938f-4141-bad7-0f34029c1eb9
using CarboKitten.Model.ALCAPS: Amount, Rate

# ╔═╡ 785153c1-2a7b-4338-bc7b-a171c4b5b5f7
using CarboKitten.Utility: in_units_of

# ╔═╡ 1bedbe94-64d5-4134-bc13-d1a0ba3056f9
const Length = typeof(1.0u"m")

# ╔═╡ e86a6fcb-5562-44e9-a40d-346fb99f226d
const Time = typeof(1.0u"Myr")

# ╔═╡ 5ff975b9-1d88-4020-b503-7e51e87de661
struct Axes
	x::Vector{Length}
	y::Vector{Length}
	bedrock_elevation::Matrix{Amount}
	t::Vector{Time}
end

# ╔═╡ b97b6042-4b64-4c13-8855-f200ab4f21a6
struct Header
	axes::Axes
	Δt::Time
	time_steps::Int
	subsidence_rate::Rate
end

# ╔═╡ 83d32f4f-49c9-45cf-b6f3-a88ac2d7c993
struct Data
	disintegration::Array{Amount, 4}
	production::Array{Amount, 4}
	deposition::Array{Amount, 4}
end

# ╔═╡ c0756361-f0a1-4e7b-8817-03a58d459ba5
function read_data(filename)
	h5open(filename) do fid
		axes = Axes(fid["input/x"][].*u"m", fid["input/y"][].*u"m", fid["input/bedrock_elevation"][].*u"m", fid["input/t"][].*u"Myr")
		attrs = HDF5.attributes(fid["input"])
		header = Header(axes, attrs["delta_t"][]*u"Myr", attrs["time_steps"][], attrs["subsidence_rate"][]*u"m/Myr")

		data = Data(fid["disintegration"][].*u"m", fid["production"][].*u"m", fid["deposition"][].*u"m")
		header, data
	end
end

# ╔═╡ b3981a26-1b82-4902-96b4-2113ecc460b3
header, data = read_data("../data/test.h5")

# ╔═╡ b051e2ce-7908-46ca-87a3-80cbabc4d426
η0 = header.axes.bedrock_elevation .- (header.axes.t[end] * header.subsidence_rate);

# ╔═╡ 02791297-3d21-470d-9de6-3132c28cba3b
η = header.axes.bedrock_elevation[:,:,[CartesianIndex()]] .+ cumsum(sum(data.deposition; dims=1)[1,:,:,:] .- sum(data.disintegration; dims=1)[1,:,:,:]; dims=3) .- (header.axes.t[end] * header.subsidence_rate);

# ╔═╡ 923d9b1f-7a79-4155-a5f9-587721718922
color = getindex.(argmax(data.deposition; dims=1)[1,:,:,:], 1);

# ╔═╡ f59b6537-f09c-4ef4-8741-61b34b1d41fe
verts = let
	verts = zeros(Float64, 100, 1001, 2)
	@views verts[:, :, 1] .= (header.axes.x |> in_units_of(u"km"))
	@views verts[:, 1, 2] = (η0[:,25] |> in_units_of(u"m"))
	@views verts[:, 2:end, 2] .= (η[:, 25, :] |> in_units_of(u"m"))
	verts
end;

# ╔═╡ 23545910-7889-429e-a21c-57897950bf3e
function explode_quad_vertices(v::Array{Float64, 3})
	w, h, _ = size(v)
	points = zeros(Float64, w, h-1, 2, 2)
	n_vertices = 2 * w * (h-1)
	n_quads = (w - 1) * (h - 1)
	@views points[:, :, 1, :] = v[1:end, 1:end-1, :]
	@views points[:, :, 2, :] = v[1:end, 2:end, :]
	idx = reshape(1:n_vertices, w, (h-1), 2)
	vtx1 = reshape(idx[1:end-1, :, 1], n_quads)
	vtx2 = reshape(idx[2:end, :, 1], n_quads)
	vtx3 = reshape(idx[2:end, :, 2], n_quads)
	vtx4 = reshape(idx[1:end-1, :, 2], n_quads)
	return reshape(points, n_vertices, 2), 
		vcat(hcat(vtx1, vtx2, vtx3), hcat(vtx1, vtx3, vtx4))
end

# ╔═╡ 0633e582-91c7-473b-815f-3fd403caca32
let
	v, f = explode_quad_vertices(verts)
	c = reshape(color[:, 25, :], 100000)
	fig = mesh(v, f, color=vcat(c, c), alpha=1.0)
	save("active-layer-test.png", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═4aaa0aee-393d-11ef-1357-739a4dfddae9
# ╠═30b991e7-5034-46c1-a46b-990f8055e4c8
# ╠═3e9d6678-1625-46c5-91c2-85f231d668a5
# ╠═3d13974e-c548-4303-bdba-c09f241fbe59
# ╠═c45863cd-1f26-4039-bd1b-37b4e7c89ec6
# ╠═59fcb379-938f-4141-bad7-0f34029c1eb9
# ╠═785153c1-2a7b-4338-bc7b-a171c4b5b5f7
# ╠═1bedbe94-64d5-4134-bc13-d1a0ba3056f9
# ╠═e86a6fcb-5562-44e9-a40d-346fb99f226d
# ╠═5ff975b9-1d88-4020-b503-7e51e87de661
# ╠═b97b6042-4b64-4c13-8855-f200ab4f21a6
# ╠═83d32f4f-49c9-45cf-b6f3-a88ac2d7c993
# ╠═c0756361-f0a1-4e7b-8817-03a58d459ba5
# ╠═b3981a26-1b82-4902-96b4-2113ecc460b3
# ╠═b051e2ce-7908-46ca-87a3-80cbabc4d426
# ╠═02791297-3d21-470d-9de6-3132c28cba3b
# ╠═923d9b1f-7a79-4155-a5f9-587721718922
# ╠═f59b6537-f09c-4ef4-8741-61b34b1d41fe
# ╠═23545910-7889-429e-a21c-57897950bf3e
# ╠═0633e582-91c7-473b-815f-3fd403caca32
