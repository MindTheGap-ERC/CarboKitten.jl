### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 4aa75138-ad54-475f-ad7e-2dfebeeb953e
begin
	using Pkg
	Pkg.activate("../../workenv")
	md"""
	!!! danger "Non-reproducible notebook"
	
		Remove this block when PR #106 is merged and a new version of CarboKitten is released.
	"""
end

# ╔═╡ 1d228208-7542-4683-8957-a0b91d5df472
using RCall

# ╔═╡ f98b493a-fdbc-11ef-1df4-f7cc750ffc21
md"""
We show how you can combine R with Julia to implement variable insolation, based on the `palinsol` R package.
"""

# ╔═╡ 759740a6-0bb4-4bdb-8e30-8600814b5752
R"install.packages('palinsol', lib='../../workenv')"

# ╔═╡ Cell order:
# ╟─4aa75138-ad54-475f-ad7e-2dfebeeb953e
# ╟─f98b493a-fdbc-11ef-1df4-f7cc750ffc21
# ╠═1d228208-7542-4683-8957-a0b91d5df472
# ╠═759740a6-0bb4-4bdb-8e30-8600814b5752
