### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ aa4723a2-8038-11ef-0a4b-a59fbd03c3ef
using Pkg; Pkg.activate("../workenv")

# ╔═╡ 19aca91a-f9af-4539-a91d-15923127cf9e
using Revise

# ╔═╡ 20390248-fe1b-4dc2-bc28-4406fd7a91f3
using ModuleMixins: @compose

# ╔═╡ 0c903628-0d63-4fcd-8e48-7351387f998b
using CarboKitten.Config: axes as box_axes

# ╔═╡ 08c22253-a895-4d5f-93b1-74f628bd6b1b
using CarboKitten.Utility: in_units_of

# ╔═╡ eb89211f-64d2-4cf7-a6ef-671dfabd4cc0
using CarboKitten.Components.Common

# ╔═╡ c9ae02f4-fee6-4742-88ae-c590fab8e03c
using CarboKitten.Model.BS92

# ╔═╡ 8c831d6f-bf47-4f4f-b2e7-12fe76776b39
using CarboKitten.Model.CAP

# ╔═╡ 89854ba4-10b9-4250-a26a-676f4012a540
using CarboKitten.Components

# ╔═╡ 47dbe7a8-78ba-4929-8ca9-14f14883fc89
using CarboKitten.Export: read_header, read_slice

# ╔═╡ ca639885-179c-4c04-9ddc-80de6d5503ec
using CarboKitten.Visualization

# ╔═╡ 741cb70c-5335-453f-95a9-d7de05ddc297
using GLMakie

# ╔═╡ deba67a0-c237-4eb3-8eb3-6594da9d03a5
using CSV

# ╔═╡ 7f91c2b5-e556-4d41-bd16-cacf5ec5545e
using DataFrames

# ╔═╡ 2a161920-6fab-4204-ad24-b921702bbe81
using Interpolations: linear_interpolation

# ╔═╡ 3adfb906-aef1-4262-89e0-ba6c318c1ebc
md"""
# Cellular Automaton
"""

# ╔═╡ e486467a-af3c-4735-b336-d989fc31f355
function run(input)
	state = CellularAutomaton.initial_state(input)
	step! = CellularAutomaton.stepper(input)
	Channel{Matrix{Int}}() do ch
		while true
			put!(ch, copy(state.ca))
			step!(state)
		end
	end
end

# ╔═╡ 54b3f849-2b3d-4c8b-b932-a4980e54c354
let
	input = CellularAutomaton.Input(
		box=Common.Box{Periodic{2}}(grid_size=(50, 50), phys_scale=1.0u"km"),
		ca_random_seed = 0,
		facies=collect(Iterators.repeated(CellularAutomaton.Facies((4, 10), (6, 10)), 3)))

	fig = Figure(size=(1000,500))
	result = run(input)
	for (i, loc) in enumerate(CartesianIndices((1:5,1:2)))
		ax = Axis(fig[loc[2], loc[1]]; title="$i")
		heatmap!(ax, take!(result))
	end
	fig
end

# ╔═╡ f6025f38-9417-43bd-ae48-f0528eb7025c
md"""
# Bosscher & Schlager
"""

# ╔═╡ c9cffcba-a08e-4769-a947-39aaf9cc1d6e
let
	function sealevel_curve()
	    data = DataFrame(CSV.File("../data/bs92-sealevel-curve.csv"))
	    linear_interpolation(data.time, data.depth)
	end
	
	input = BS92.Input(
	    tag = "example model BS92",
	    box = Common.Box{Shelf}(grid_size=(100, 1), phys_scale=600.0u"m"),
	    time = TimeProperties(
	      Δt = 10.0u"yr",
	      steps = 8000,
	      write_interval = 100),
	    sea_level = let sc = sealevel_curve()
	      t -> -sc(t / u"yr") * u"m"
	    end,
	    bedrock_elevation = (x, y) -> - x / 300.0,
	    subsidence_rate = 0.0u"m/yr",
	    insolation = 400.0u"W/m^2",
	    facies = [BS92.Facies(
	      maximum_growth_rate = 0.005u"m/yr",
	      saturation_intensity = 50.0u"W/m^2",
	      extinction_coefficient = 0.05u"m^-1"
	    )])

	H5Writer.run(Model{BS92}, input, "../data/output/bs92.h5")
end

# ╔═╡ f23fa7c3-8bdd-43af-af1a-f6fb36bdd92f
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

# ╔═╡ e93b7ac8-dcf9-4ca1-8d1b-01d6dd139930
md"""
# Multiple facies
"""

# ╔═╡ bac48b9d-46d0-44a4-9e9f-7fa32b8fc75d
facies = [
	BS92.Facies(
        maximum_growth_rate=500u"m/Myr"/4,
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2"),
    BS92.Facies(
        maximum_growth_rate=400u"m/Myr"/4,
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2"),
    BS92.Facies(
        maximum_growth_rate=100u"m/Myr"/4,
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2")
]

# ╔═╡ 79b3c99b-0691-405c-8ed8-1da6d25cb9a8
let
	input = BS92.Input(
	    tag = "example model BS92",
	    box = Common.Box{Shelf}(grid_size=(100, 1), phys_scale=150.0u"m"),
	    time = TimeProperties(
	      Δt = 200.0u"yr",
	      steps = 5000,
	      write_interval = 1),
	    sea_level = t -> 4.0u"m" * sin(2π * t / 0.2u"Myr"),
	    bedrock_elevation = (x, y) -> - x / 300.0,
	    subsidence_rate = 50.0u"m/Myr",
	    insolation = 400.0u"W/m^2",
	    facies = facies)

	H5Writer.run(Model{BS92}, input, "../data/output/bs92-3facies.h5")
end

# ╔═╡ 67f405a6-1252-4fc8-9c70-0d97dd4437ae
let
	header, data = read_slice("../data/output/bs92-3facies.h5", :, 1)
	sediment_profile(header, data)
end

# ╔═╡ 20901dd1-ce07-4dd0-aa1f-3673a8012da3
md"""
# Production + CA
"""

# ╔═╡ a864d6ea-99bf-4513-9f44-f15b1d19d6af
let
	cap_facies = [
	    CAP.Facies(
			viability_range = (4, 10),
			activation_range = (6, 10),
			maximum_growth_rate = 500u"m/Myr",
			extinction_coefficient = 0.8u"m^-1",
			saturation_intensity = 60u"W/m^2"),
	
	    CAP.Facies(
			viability_range = (4, 10),
			activation_range = (6, 10),
			maximum_growth_rate = 400u"m/Myr",
			extinction_coefficient = 0.1u"m^-1",
			saturation_intensity = 60u"W/m^2"),
	
	    CAP.Facies(
			viability_range = (4, 10),
			activation_range = (6, 10),
			maximum_growth_rate = 100u"m/Myr",
			extinction_coefficient = 0.005u"m^-1",
			saturation_intensity = 60u"W/m^2")
	]

	input = CAP.Input(
		tag = "example model BS92",
		box = Common.Box{Shelf}(grid_size=(100, 50), phys_scale=150.0u"m"),
		time = TimeProperties(
			Δt = 200.0u"yr",
			steps = 5000,
			write_interval = 10),
		sea_level = t -> 4.0u"m" * sin(2π * t / 0.2u"Myr"),
		bedrock_elevation = (x, y) -> - x / 300.0,
		subsidence_rate = 50.0u"m/Myr",
		insolation = 400.0u"W/m^2",
		facies = cap_facies)

	H5Writer.run(Model{CAP}, input, "../data/output/cap1.h5")
end

# ╔═╡ b523780a-7841-4f89-b144-f7b277bec831
let
	header, data = read_slice("../data/output/cap1.h5", :, 25)
	sediment_profile(header, data)
end

# ╔═╡ Cell order:
# ╠═aa4723a2-8038-11ef-0a4b-a59fbd03c3ef
# ╠═19aca91a-f9af-4539-a91d-15923127cf9e
# ╠═20390248-fe1b-4dc2-bc28-4406fd7a91f3
# ╠═0c903628-0d63-4fcd-8e48-7351387f998b
# ╠═08c22253-a895-4d5f-93b1-74f628bd6b1b
# ╠═eb89211f-64d2-4cf7-a6ef-671dfabd4cc0
# ╠═c9ae02f4-fee6-4742-88ae-c590fab8e03c
# ╠═8c831d6f-bf47-4f4f-b2e7-12fe76776b39
# ╠═89854ba4-10b9-4250-a26a-676f4012a540
# ╠═47dbe7a8-78ba-4929-8ca9-14f14883fc89
# ╠═ca639885-179c-4c04-9ddc-80de6d5503ec
# ╠═741cb70c-5335-453f-95a9-d7de05ddc297
# ╠═deba67a0-c237-4eb3-8eb3-6594da9d03a5
# ╠═7f91c2b5-e556-4d41-bd16-cacf5ec5545e
# ╠═2a161920-6fab-4204-ad24-b921702bbe81
# ╟─3adfb906-aef1-4262-89e0-ba6c318c1ebc
# ╠═e486467a-af3c-4735-b336-d989fc31f355
# ╠═54b3f849-2b3d-4c8b-b932-a4980e54c354
# ╟─f6025f38-9417-43bd-ae48-f0528eb7025c
# ╠═c9cffcba-a08e-4769-a947-39aaf9cc1d6e
# ╠═f23fa7c3-8bdd-43af-af1a-f6fb36bdd92f
# ╟─e93b7ac8-dcf9-4ca1-8d1b-01d6dd139930
# ╠═bac48b9d-46d0-44a4-9e9f-7fa32b8fc75d
# ╠═79b3c99b-0691-405c-8ed8-1da6d25cb9a8
# ╠═67f405a6-1252-4fc8-9c70-0d97dd4437ae
# ╟─20901dd1-ce07-4dd0-aa1f-3673a8012da3
# ╠═a864d6ea-99bf-4513-9f44-f15b1d19d6af
# ╠═b523780a-7841-4f89-b144-f7b277bec831
