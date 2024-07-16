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

# ╔═╡ 6be3c53b-88a4-4634-b7a4-050c99824d15
using Statistics: mean

# ╔═╡ 1bedbe94-64d5-4134-bc13-d1a0ba3056f9
const Length = typeof(1.0u"m")

# ╔═╡ e86a6fcb-5562-44e9-a40d-346fb99f226d
const Time = typeof(1.0u"Myr")

# ╔═╡ 8658fcc2-1545-45c3-85a0-ad0ec543e58d
na = [CartesianIndex()]

# ╔═╡ 5ff975b9-1d88-4020-b503-7e51e87de661
struct Axes
	x::Vector{Length}
	y::Vector{Length}
	t::Vector{Time}
end

# ╔═╡ b97b6042-4b64-4c13-8855-f200ab4f21a6
struct Header
	axes::Axes
	Δt::Time
	time_steps::Int
	bedrock_elevation::Matrix{Amount}
	sea_level::Vector{Length}
	subsidence_rate::Rate
end

# ╔═╡ 83d32f4f-49c9-45cf-b6f3-a88ac2d7c993
struct Data
	disintegration::Array{Amount, 4}
	production::Array{Amount, 4}
	deposition::Array{Amount, 4}
	sediment_elevation::Array{Amount, 3}
end

# ╔═╡ 1a439db1-3d14-4678-82cc-63f2b975abd7
struct DataSlice
	disintegration::Array{Amount, 3}
	production::Array{Amount, 3}
	deposition::Array{Amount, 3}
	sediment_elevation::Array{Amount, 2}
end

# ╔═╡ 775bc869-2d92-4850-89dd-2af588d44ae6
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

# ╔═╡ c0756361-f0a1-4e7b-8817-03a58d459ba5
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

# ╔═╡ 458a4d96-fa19-48c6-bf9c-a6c01760f29d
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

# ╔═╡ ba7d41ba-20fe-4ce6-92d2-625b893be1f7
read_slice("../data/test.h5", :, :, 25, :)

# ╔═╡ b3981a26-1b82-4902-96b4-2113ecc460b3
header, data = read_data("../data/test.h5")

# ╔═╡ b051e2ce-7908-46ca-87a3-80cbabc4d426
η0 = header.bedrock_elevation .- (header.axes.t[end] * header.subsidence_rate);

# ╔═╡ 6f803392-8a12-41f9-9085-35ce9555323e
elevation(h::Header, d::Data) = let bl = h.bedrock_elevation[:,:,na],
									sr = h.axes.t[end] * h.subsidence_rate
	cat(bl, bl .+ d.sediment_elevation; dims=3) .- sr
end

# ╔═╡ e8307d25-1432-40cb-b89e-e1005631cecb
elevation(h::Header, d::DataSlice, y) = let bl = h.bedrock_elevation[:,y,na],
	       								 sr = h.axes.t[end] * h.subsidence_rate
	cat(bl, bl .+ d.sediment_elevation; dims=2) .- sr
end

# ╔═╡ 5bca1081-fb71-475a-be45-1551ee05861c
η = elevation(header, data);

# ╔═╡ 923d9b1f-7a79-4155-a5f9-587721718922
colormax(d::Data) = getindex.(argmax(d.deposition; dims=1)[1,:,:,:], 1)

# ╔═╡ 3e9f7a26-48be-4f86-8c6f-bbdce56b84ab
colormax(d::DataSlice) = getindex.(argmax(d.deposition; dims=1)[1,:,:], 1)

# ╔═╡ d2f725e6-dfb5-43a6-9328-ba3662205757
colorrgb = eachslice(replace(data.deposition ./ sum(data.deposition; dims=1), NaN=>0.0), dims=(2,3,4)) .|> splat(RGBf);

# ╔═╡ 23545910-7889-429e-a21c-57897950bf3e
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

# ╔═╡ 2f9068db-7feb-4652-ba87-ff0facb88e68
md"""
## Production profile
"""

# ╔═╡ 669ba75e-803a-4ea7-b75c-6e375f30321b
let
	y = 25
	prod = sum(data.production; dims=1)[1,:,:,:]
	colormax = getindex.(argmax(data.production; dims=1)[1,:,:,:], 1)
	c = Array{Union{Int, Missing}}(missing, size(prod)...)
	mask = prod .!= 0.0u"m"
	c[mask] .= colormax[mask]
	heatmap(header.axes.x |> in_units_of(u"km"), header.axes.t[1:end-1] |> in_units_of(u"Myr"), c[:,y,:], colormap=:viridis)
end

# ╔═╡ 46bd4f32-4fba-403d-8d2e-58397379f06d
md"""
## Detect production gaps
"""

# ╔═╡ 2603c0fb-a333-48a0-ac8c-de04c37bb867
local_height = η .- (header.subsidence_rate .* (header.axes.t .- header.axes.t[end]) .+ header.sea_level)[na, na, :];

# ╔═╡ 465ae61e-9b18-41ed-9aca-acf46a15fd67
let
	fig = Figure()
	ax = Axis(fig[1, 1])
	heatmap!(ax, local_height[:,25,:] .> 0u"m")
	fig
end

# ╔═╡ 5bf17772-37bf-4272-9648-ef716c6bd89d
let
	x = header.axes.x |> in_units_of(u"km")
	t = header.axes.t |> in_units_of(u"Myr")
	lh = η .- (header.subsidence_rate .* (header.axes.t .- header.axes.t[end]) .+ header.sea_level)[na,na, :]
	fig = Figure(size=(800, 400))
	ax = Axis3(fig[1,1:3], azimuth=7π/4, zlabel="", xlabel="x(km)", ylabel="t(Myr)")
	surface!(ax, x, t, lh[:,25,:] |> in_units_of(u"m"))
	ax2 = Axis(fig[1,4:5], xlabel="t(Myr)", ylabel="-w(m)", title="x=1km")
	lines!(ax2, t, lh[8,25,:] |> in_units_of(u"m"))
	fig
end

# ╔═╡ 69f4ccf2-e204-4bc8-906e-d69bc927760e
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

# ╔═╡ 42787e9b-9037-487c-b551-609fa0e343f7
bean_counter(BitArray([0 1 0; 0 0 0; 1 1 0]))

# ╔═╡ a41be475-5632-4d62-89ce-7810efeaab12
gaps, n_gaps = bean_counter(local_height .> 0u"m")

# ╔═╡ d7abbfb1-020b-4a09-8df6-7adf698e7388
heatmap(gaps[:, 25, :] .|> (x-> x==0 ? missing : x))

# ╔═╡ 66c0f8c8-6336-45ea-9fca-db7a4638efd6
md"""
## Sediment profile
"""

# ╔═╡ 0633e582-91c7-473b-815f-3fd403caca32
function plot_sediment_profile(filename, y)
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

	fig = Figure(size=(800,600))
	ax = Axis(fig[1,1])
	c = reshape(colormax(data)[:, :], length(x) * (length(t) - 1))
	mesh!(ax, v, f, color=vcat(c, c), alpha=1.0)

	for g = 1:n_gaps
		size = sum(gaps .== g)
		if size < 1000
			continue
		end
		gap = mean.(skipmissing.(eachslice(CartesianIndices(ξ) .|> (i -> gaps[i] == g ? ξ[i] : missing), dims=(1,))))
		lines!(ax, x, gap |> in_units_of(u"m"), color=:white, linewidth=2, linestyle=:dash)
	end
	
	fig
end

# ╔═╡ 79fdebf7-740c-4bce-945d-38c997ce3b9b
plot_sediment_profile("../data/test.h5", 25)

# ╔═╡ 70affadb-9eed-4bb1-9515-b4625a75ff70
md"""
## Sediment profile with time

Seen from above this is a Wheeler diagram, from the side we get a sediment profile.
"""

# ╔═╡ 2029a364-3ad1-47f8-ace3-e3698f1f557e
let
	y = 25
	x = header.axes.x |> in_units_of(u"km")
	t = header.axes.t |> in_units_of(u"Myr")

	verts = zeros(Float64, length(x), length(t), 3)
	@views verts[:, :, 1] .= (header.axes.x |> in_units_of(u"km"))
	@views verts[:, :, 3] .= (η[:, y, :] |> in_units_of(u"m"))
	@views verts[:, :, 2] .= (header.axes.t |> in_units_of(u"Myr"))[na,:]
	v, f = explode_quad_vertices(verts)
	c = reshape(colormax(data)[:, y, :], length(x) * (length(t) - 1))

	fig = Figure()
	ax = Axis3(fig[1,1], azimuth=7π/4, xlabel="x(km)", ylabel="t(Myr)", zlabel="z(m)")
	mesh!(ax, v, f, color=vcat(c, c), alpha=1.0)
	save("active-layer-test.png", fig)
	fig
end

# ╔═╡ e151095a-78c8-4e12-af2e-f8c85a6cf51f
md"""
## Top layer
"""

# ╔═╡ 9651674a-bd63-49dd-ad34-5357ee1d6caf
let
	fig = Figure()
	ax = Axis3(fig[1, 1], azimuth=7π/4)
	surface!(ax, header.axes.x |> in_units_of(u"km"), header.axes.y |> in_units_of(u"km"), η[:,:,end] |> in_units_of(u"m"), color=colorrgb[:,:,end])
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
# ╠═8658fcc2-1545-45c3-85a0-ad0ec543e58d
# ╠═5ff975b9-1d88-4020-b503-7e51e87de661
# ╠═b97b6042-4b64-4c13-8855-f200ab4f21a6
# ╠═83d32f4f-49c9-45cf-b6f3-a88ac2d7c993
# ╠═1a439db1-3d14-4678-82cc-63f2b975abd7
# ╠═775bc869-2d92-4850-89dd-2af588d44ae6
# ╠═c0756361-f0a1-4e7b-8817-03a58d459ba5
# ╠═458a4d96-fa19-48c6-bf9c-a6c01760f29d
# ╠═ba7d41ba-20fe-4ce6-92d2-625b893be1f7
# ╠═b3981a26-1b82-4902-96b4-2113ecc460b3
# ╠═b051e2ce-7908-46ca-87a3-80cbabc4d426
# ╠═6f803392-8a12-41f9-9085-35ce9555323e
# ╠═5bca1081-fb71-475a-be45-1551ee05861c
# ╠═e8307d25-1432-40cb-b89e-e1005631cecb
# ╠═923d9b1f-7a79-4155-a5f9-587721718922
# ╠═3e9f7a26-48be-4f86-8c6f-bbdce56b84ab
# ╠═d2f725e6-dfb5-43a6-9328-ba3662205757
# ╠═23545910-7889-429e-a21c-57897950bf3e
# ╟─2f9068db-7feb-4652-ba87-ff0facb88e68
# ╠═669ba75e-803a-4ea7-b75c-6e375f30321b
# ╟─46bd4f32-4fba-403d-8d2e-58397379f06d
# ╠═6be3c53b-88a4-4634-b7a4-050c99824d15
# ╠═2603c0fb-a333-48a0-ac8c-de04c37bb867
# ╠═465ae61e-9b18-41ed-9aca-acf46a15fd67
# ╠═5bf17772-37bf-4272-9648-ef716c6bd89d
# ╠═69f4ccf2-e204-4bc8-906e-d69bc927760e
# ╠═42787e9b-9037-487c-b551-609fa0e343f7
# ╠═a41be475-5632-4d62-89ce-7810efeaab12
# ╠═d7abbfb1-020b-4a09-8df6-7adf698e7388
# ╟─66c0f8c8-6336-45ea-9fca-db7a4638efd6
# ╠═0633e582-91c7-473b-815f-3fd403caca32
# ╠═79fdebf7-740c-4bce-945d-38c997ce3b9b
# ╟─70affadb-9eed-4bb1-9515-b4625a75ff70
# ╠═2029a364-3ad1-47f8-ace3-e3698f1f557e
# ╟─e151095a-78c8-4e12-af2e-f8c85a6cf51f
# ╠═9651674a-bd63-49dd-ad34-5357ee1d6caf
