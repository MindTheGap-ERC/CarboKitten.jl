### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 74ccac56-7b4d-11ef-3d4a-01e35a06c040
using Pkg; Pkg.activate("../workenv")

# ╔═╡ d9671eec-2705-48cd-b019-1a1a34cd8cb4
using Revise

# ╔═╡ 4386d059-ad7b-4b61-801d-85a80055549b
using CairoMakie

# ╔═╡ aa642ec8-9ab3-4a2d-851f-6295538351c0
using CarboKitten.Export: read_slice, Header, DataSlice

# ╔═╡ 4f5ae7d3-bfc1-429a-8f6e-875cbf32d784
using CarboKitten.Stencil: convolution

# ╔═╡ d79f908b-bace-4a11-95ea-c3a4d70fb5f9
using CarboKitten.BoundaryTrait

# ╔═╡ 935346d3-e8a6-4b6b-97f5-70cbf07964d2
using Unitful

# ╔═╡ efff0b3f-c3c0-4e7b-a46a-840227646095
using CarboKitten.Visualization: sediment_profile, wheeler_diagram

# ╔═╡ 3b39f5f4-c7ae-47c0-ab67-79484d38af67
header, data = read_slice("../data/alcaps2.h5", :, 25)

# ╔═╡ 4fc3ec1b-596e-4bff-8717-ce119a10e0ba
const na = [CartesianIndex()]

# ╔═╡ 1d287bac-e53c-4605-b4c7-e52c031a338b
elevation(h::Header, d::DataSlice) =
    let bl = h.bedrock_elevation[d.slice..., na],
        sr = h.axes.t[end] * h.subsidence_rate

        cat(bl, bl .+ d.sediment_elevation; dims=2) .- sr
    end

# ╔═╡ 2a12beaf-5139-4771-98e5-bb3eb3cd5768
water_depth(header::Header, data::DataSlice) =
    let h = elevation(header, data),
        s = header.subsidence_rate .* (header.axes.t .- header.axes.t[end]),
        l = header.sea_level

        h .- (s .+ header.sea_level)[na, :]
    end

# ╔═╡ c17123c4-8998-4805-a248-e51a4d197efd
colormax(d) = getindex.(argmax(d; dims=1)[1, :, :], 1)

# ╔═╡ 18642e2a-36ed-44f1-9f24-16d3a7428433
let
	fig = sediment_profile(header, data)
	save("sediment_profile.png", fig)
	fig
end

# ╔═╡ cc1f98a9-b689-4ce0-a1f1-9d0263ebec63
const Depth = typeof(1.0u"m")

# ╔═╡ 01839fae-7756-4be4-8017-c32aa03dba0d
let
	fig = wheeler_diagram(header, data; smooth_size=(3,11))
	save("wheeler_test_noblur.png", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═74ccac56-7b4d-11ef-3d4a-01e35a06c040
# ╠═d9671eec-2705-48cd-b019-1a1a34cd8cb4
# ╠═4386d059-ad7b-4b61-801d-85a80055549b
# ╠═aa642ec8-9ab3-4a2d-851f-6295538351c0
# ╠═4f5ae7d3-bfc1-429a-8f6e-875cbf32d784
# ╠═d79f908b-bace-4a11-95ea-c3a4d70fb5f9
# ╠═935346d3-e8a6-4b6b-97f5-70cbf07964d2
# ╠═efff0b3f-c3c0-4e7b-a46a-840227646095
# ╠═3b39f5f4-c7ae-47c0-ab67-79484d38af67
# ╠═4fc3ec1b-596e-4bff-8717-ce119a10e0ba
# ╠═1d287bac-e53c-4605-b4c7-e52c031a338b
# ╠═2a12beaf-5139-4771-98e5-bb3eb3cd5768
# ╠═c17123c4-8998-4805-a248-e51a4d197efd
# ╠═18642e2a-36ed-44f1-9f24-16d3a7428433
# ╠═cc1f98a9-b689-4ce0-a1f1-9d0263ebec63
# ╠═01839fae-7756-4be4-8017-c32aa03dba0d
