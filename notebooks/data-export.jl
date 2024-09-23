### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ e84f6252-7043-11ef-19ec-cbda366d00ef
using Pkg; Pkg.activate("../workenv")

# ╔═╡ ef15eb8e-e70a-4c0c-9037-eae7cd65515d
using Revise

# ╔═╡ 32dd37cc-e4db-43d7-90d0-e197e552fda3
using CarboKitten.Export: Header, Data, data_export, read_data, extract_sac, CSV, age_depth_model, extract_sc

# ╔═╡ ae90ebbd-2658-400b-89c2-9acaa1d676a9
using DataFrames

# ╔═╡ a4b39575-11c9-404d-9333-3296c6052902
using Unitful

# ╔═╡ 5201a517-aa6c-482c-af02-7bef3d2a903f
using GLMakie

# ╔═╡ 1780ccbd-9c49-455e-8987-a4cacb0d05fc
using CarboKitten.Visualization: sediment_profile!

# ╔═╡ 777f6eea-3901-492f-8fa2-0e9508594b1a
header, data = read_data("../data/alcaps_default.h5")

# ╔═╡ b6ec79cd-0d67-48f0-9111-58967098f7d2
open("../data/export_test_sc.csv", "r") do io
	print(read(io, String))
end

# ╔═╡ 51f80e72-a691-4132-b8dd-9e7ac8a2f165
grid_locations=[(10, 25), (30, 25), (50, 25), (70, 25)]

# ╔═╡ 8d67a359-69b2-4d08-a8bc-a17b2c8975af
data_export(CSV(grid_locations=grid_locations,
	sac = "../data/export_test_sac.csv",
	adm = "../data/export_test_adm.csv",
	sc = "../data/export_test_sc.csv"), header, data)

# ╔═╡ 58776042-49a1-4d31-af06-366066b9c019
sac = extract_sac(header, data, grid_locations)

# ╔═╡ 2189c89d-cc17-4fe6-94d0-af16500a6d61
replace.(filter(contains("sac"), names(sac)), "sac" => "adm")

# ╔═╡ 328d1f1b-8865-49e6-abd3-01dbebe3c0c0
age_depth_model(sac)

# ╔═╡ 8bd0e45d-3c9e-4e38-8ec7-3d70d2f5ff0d
let
	fig = Figure()
	ax = Axis(fig[1,1])
	sediment_profile!(ax, "../data/alcaps_default.h5", 25)
	vlines!(ax, header.axes.x[getindex.(grid_locations, 1)] / u"km")
	fig
end

# ╔═╡ 6a1ed82a-d00f-4bdd-b8d3-899efbceabf6
let
	fig = Figure()
	ax = Axis(fig[1,1])
	for (i, g) in enumerate(grid_locations)
		lines!(ax, header.axes.t[1:end-1] / u"Myr", sac[!,"sac$i"] / u"m")
		lines!(ax, header.axes.t[1:end-1] / u"Myr", age_depth_model(sac[!,"sac$i"]) / u"m")
	end
	fig
end

# ╔═╡ 0df15e29-e66a-4904-808d-0442620d6d00
columns(i, loc) = let n_facies = size(data.production)[1]
	("sac$(i)" => data.sediment_elevation[loc..., :],
	Iterators.flatten(("fac$(f)_dis$(i)" => data.disintegration[f, loc..., :],
	  "fac$(f)_prod$(i)" => data.production[f, loc..., :],
	  "fac$(f)_dep$(f)" => data.deposition[f, loc..., :])
	 for f in 1:n_facies)...)
end


# ╔═╡ 49ec1a8f-e6a5-4084-9e42-7448e3d1fcc1
function stratigraphic_column(header::Header, data::Data, loc::NTuple{2,Int})
	n_facies = size(data.production)[1]
	n_times = length(header.axes.t) - 1 
	sc = zeros(typeof(1.0u"m"), n_facies, n_times)

	for f = 1:n_facies
	for ts = 1:n_times
		acc = data.deposition[f, loc..., ts] - data.disintegration[f, loc..., ts]
		if acc > 0.0u"m"
			sc[f, ts] = acc
			continue
		end
		ts_down = ts - 1
		while acc < 0.0u"m"
			ts_down < 1 && break
			if -acc < sc[f, ts_down]
				sc[f, ts_down] -= acc
				break
			else
				acc += sc[f, ts_down]
				sc[f, ts_down] = 0.0u"m"
			end
			ts_down -= 1
		end
	end
	end

	sc
end

# ╔═╡ 19fa2885-3886-425e-8e34-01ffa40518ad
extract_sc(header, data, grid_locations)

# ╔═╡ 1190809d-a887-4f05-961f-fc20c3434ce9
[1, 2, 3] .=> [:a, :b, :c]

# ╔═╡ b8968897-9a60-4b21-9ca0-ec382bc30917
DataFrame(stratigraphic_column(header, data, (40, 25))', ["sc1", "sc2", "sc3"])

# ╔═╡ b1fa959e-a364-40ca-811a-53ea29543eed
let
	sc = stratigraphic_column(header, data, grid_locations[1])
	fig = Figure()
	ax = Axis(fig[1, 1])
	for i = 1:3
		lines!(ax, header.axes.t[1:end-1] / u"Myr", cumsum(sc[i,:]) / u"m")
	end
	lines!(ax, header.axes.t[1:end-1] / u"Myr", sac[!,2] / u"m")
	fig
end

# ╔═╡ Cell order:
# ╠═e84f6252-7043-11ef-19ec-cbda366d00ef
# ╠═ef15eb8e-e70a-4c0c-9037-eae7cd65515d
# ╠═32dd37cc-e4db-43d7-90d0-e197e552fda3
# ╠═ae90ebbd-2658-400b-89c2-9acaa1d676a9
# ╠═a4b39575-11c9-404d-9333-3296c6052902
# ╠═5201a517-aa6c-482c-af02-7bef3d2a903f
# ╠═1780ccbd-9c49-455e-8987-a4cacb0d05fc
# ╠═777f6eea-3901-492f-8fa2-0e9508594b1a
# ╠═8d67a359-69b2-4d08-a8bc-a17b2c8975af
# ╠═b6ec79cd-0d67-48f0-9111-58967098f7d2
# ╠═51f80e72-a691-4132-b8dd-9e7ac8a2f165
# ╠═58776042-49a1-4d31-af06-366066b9c019
# ╠═2189c89d-cc17-4fe6-94d0-af16500a6d61
# ╠═328d1f1b-8865-49e6-abd3-01dbebe3c0c0
# ╠═8bd0e45d-3c9e-4e38-8ec7-3d70d2f5ff0d
# ╠═6a1ed82a-d00f-4bdd-b8d3-899efbceabf6
# ╠═0df15e29-e66a-4904-808d-0442620d6d00
# ╠═49ec1a8f-e6a5-4084-9e42-7448e3d1fcc1
# ╠═19fa2885-3886-425e-8e34-01ffa40518ad
# ╠═1190809d-a887-4f05-961f-fc20c3434ce9
# ╠═b8968897-9a60-4b21-9ca0-ec382bc30917
# ╠═b1fa959e-a364-40ca-811a-53ea29543eed
