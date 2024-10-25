### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 38e3126f-441f-475a-b1aa-0485dd1831bc
using Pkg; Pkg.activate("../workenv")

# ╔═╡ 1d7ba24e-8bcc-11ef-00f2-35768ab3077f
using DelimitedFiles

# ╔═╡ 423c80c7-fcf3-4f21-b474-32b8da52bc44
using DataFrames

# ╔═╡ 3574830d-bbaf-4f5f-a80a-c8b2a88cf93b
using Interpolations

# ╔═╡ aae1e4d0-8035-464f-9257-e312d70d9976
using GLMakie

# ╔═╡ 070b6d00-f74d-471f-899e-ef27b935b7e4
using Unitful

# ╔═╡ ec2481a6-92e4-4b33-ba46-72936117ee30


# ╔═╡ 7a1e0de4-f88e-4720-8540-2c269d7a8a47
data, header = readdlm("../data/Cenozoic_sea_level_reconstruction.tab", '\t', skipstart=40, header=true)

# ╔═╡ 06a39be4-3200-49dc-af95-6ae6fd640462
header

# ╔═╡ 7e2c9e35-b922-422c-b7ef-6b4d2bc1bd6f
df = DataFrame(age=data[:, 4]*u"kyr", sealevel=data[:, 7]*u"m")

# ╔═╡ bed095c2-636e-4fe4-aa0e-43275c96f805
let
	fig = Figure()
	ax = Axis(fig[1,1], xreversed=true,
		dim1_conversion=Makie.UnitfulConversion(u"Myr"),
		limits=((0, 1.0), nothing))
	lines!(ax, df.age, df.sealevel)
	fig
end

# ╔═╡ 2168cd4a-8707-41e5-b017-d33334e04322
linear_interpolation(df.age, df.sealevel)

# ╔═╡ Cell order:
# ╠═38e3126f-441f-475a-b1aa-0485dd1831bc
# ╠═1d7ba24e-8bcc-11ef-00f2-35768ab3077f
# ╠═423c80c7-fcf3-4f21-b474-32b8da52bc44
# ╠═3574830d-bbaf-4f5f-a80a-c8b2a88cf93b
# ╠═aae1e4d0-8035-464f-9257-e312d70d9976
# ╠═070b6d00-f74d-471f-899e-ef27b935b7e4
# ╠═ec2481a6-92e4-4b33-ba46-72936117ee30
# ╠═7a1e0de4-f88e-4720-8540-2c269d7a8a47
# ╠═06a39be4-3200-49dc-af95-6ae6fd640462
# ╠═7e2c9e35-b922-422c-b7ef-6b4d2bc1bd6f
# ╠═bed095c2-636e-4fe4-aa0e-43275c96f805
# ╠═2168cd4a-8707-41e5-b017-d33334e04322
