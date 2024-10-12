### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 8433507c-88b1-11ef-201c-df0d6533a9fe
using Pkg; Pkg.activate("../workenv")

# ╔═╡ af8f7212-a873-4d76-8a7a-e876c1ebb6e5
using Revise

# ╔═╡ e671c760-0c57-432b-8cc1-4c7e85b8b514
using GLMakie

# ╔═╡ 9ddaa447-5f7c-4293-b5bc-668857d8ad16
using CarboKitten.Visualization

# ╔═╡ c9be5be5-e22b-498c-af4b-3b901a6edac6
using CarboKitten.Export: read_header

# ╔═╡ 5411d93f-2c76-40eb-998a-fa315090f45f
using HDF5

# ╔═╡ cee9e8a5-d897-42ba-a2b1-53bee950e366
summary_plot("../data/output/alcap2.h5"; wheeler_smooth=(3, 11))

# ╔═╡ 5fb09723-5d14-4132-ba68-809a27fade6e
summary_plot("../data/output/bs92.h5")

# ╔═╡ 541a3f3f-d02b-4382-ad75-4f27872114c4
summary_plot("../data/output/cap1.h5")

# ╔═╡ 28a09fe2-edb6-4d3e-b20a-f8a86efcbd38
h5open("../data/output/bs92.h5", "r") do fid
	h = read_header(fid)
	minimum(h.bedrock_elevation)
end

# ╔═╡ Cell order:
# ╠═8433507c-88b1-11ef-201c-df0d6533a9fe
# ╠═af8f7212-a873-4d76-8a7a-e876c1ebb6e5
# ╠═e671c760-0c57-432b-8cc1-4c7e85b8b514
# ╠═9ddaa447-5f7c-4293-b5bc-668857d8ad16
# ╠═c9be5be5-e22b-498c-af4b-3b901a6edac6
# ╠═5411d93f-2c76-40eb-998a-fa315090f45f
# ╠═cee9e8a5-d897-42ba-a2b1-53bee950e366
# ╠═5fb09723-5d14-4132-ba68-809a27fade6e
# ╠═541a3f3f-d02b-4382-ad75-4f27872114c4
# ╠═28a09fe2-edb6-4d3e-b20a-f8a86efcbd38