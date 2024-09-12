### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ e84f6252-7043-11ef-19ec-cbda366d00ef
using Pkg; Pkg.activate("../workenv")

# ╔═╡ ef15eb8e-e70a-4c0c-9037-eae7cd65515d
using Revise

# ╔═╡ 32dd37cc-e4db-43d7-90d0-e197e552fda3
using CarboKitten.Export: data_export, read_data, extract_sac

# ╔═╡ ae90ebbd-2658-400b-89c2-9acaa1d676a9
using DataFrames

# ╔═╡ 4efdc2d1-8dd3-421a-9792-977911033d2b
using CSV

# ╔═╡ a4b39575-11c9-404d-9333-3296c6052902
using Unitful

# ╔═╡ 33de7b7e-e13d-432d-9b58-413d289528ef


# ╔═╡ 777f6eea-3901-492f-8fa2-0e9508594b1a
header, data = read_data("../data/alcaps_default.h5")

# ╔═╡ 8d67a359-69b2-4d08-a8bc-a17b2c8975af
df = extract_sac(header, data, [(10, 25)])

# ╔═╡ 0f953d8b-d18d-47ce-b5f6-3dbddfa27da9
function unitful_headers(df::DataFrame)
	["$(e.variable) [$(unit(e.eltype))]" for e in eachrow(describe(df))]
end

# ╔═╡ 100ea223-4af5-42ea-8ede-9a18579862e8
function Unitful.ustrip(df::DataFrame)
	DataFrame((e.variable => df[!,e.variable] / unit(e.eltype)
		for e in eachrow(describe(df)))...)
end

# ╔═╡ 08c6dc90-f1a3-48bb-899e-8920ce1e241d
function write_unitful_csv(io::IO, df::DataFrame)
	CSV.write(io, ustrip(df), header=unitful_headers(df))
end

# ╔═╡ 73a032a1-8626-4abd-8eab-42814859b41d
let
	buffer = UInt8[]
	io = IOBuffer(buffer, write=true)
	write_unitful_csv(io, df)
	print(read(IOBuffer(buffer), String))
end

# ╔═╡ Cell order:
# ╠═e84f6252-7043-11ef-19ec-cbda366d00ef
# ╠═ef15eb8e-e70a-4c0c-9037-eae7cd65515d
# ╠═32dd37cc-e4db-43d7-90d0-e197e552fda3
# ╠═ae90ebbd-2658-400b-89c2-9acaa1d676a9
# ╠═4efdc2d1-8dd3-421a-9792-977911033d2b
# ╠═a4b39575-11c9-404d-9333-3296c6052902
# ╠═33de7b7e-e13d-432d-9b58-413d289528ef
# ╠═777f6eea-3901-492f-8fa2-0e9508594b1a
# ╠═8d67a359-69b2-4d08-a8bc-a17b2c8975af
# ╠═0f953d8b-d18d-47ce-b5f6-3dbddfa27da9
# ╠═100ea223-4af5-42ea-8ede-9a18579862e8
# ╠═08c6dc90-f1a3-48bb-899e-8920ce1e241d
# ╠═73a032a1-8626-4abd-8eab-42814859b41d
