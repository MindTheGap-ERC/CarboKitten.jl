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

# ╔═╡ bc3cc455-f7cb-4c57-9e7a-f78ae55afe14
using CarboKitten.Model: ALCAP2, CAP, BS92

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

# ╔═╡ ae356c01-2520-4bdf-8555-8b465040e494
let
	header, data = read_slice("../data/output/bs92.h5", :, 1)
	sediment_profile(header, data)
end

# ╔═╡ c954e223-fb91-4f86-8199-f91f39156773
let
	header, data = read_slice("../data/output/bs92.h5", :, 1)
	wheeler_diagram(header, data; smooth_size=(1,1), range=(-5000.0u"m/Myr", 5000.0u"m/Myr"))
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

# ╔═╡ c8fbbf2a-de9e-4e88-b3bd-201091e8d869


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

# ╔═╡ 72e1d116-06f0-4623-b068-9dd32c6f0eba
md"""
# Active Layer Transport
"""

# ╔═╡ c3b7220f-49f4-42a2-bd2a-1c554e86d3c8
let
	facies = [
    	ALCAP2.Facies(
			viability_range = (4, 10),
            activation_range = (6, 10),
            maximum_growth_rate = 500u"m/Myr",
            extinction_coefficient = 0.8u"m^-1",
            saturation_intensity = 60u"W/m^2",
            diffusion_coefficient = 1000u"m"),

     	ALCAP2.Facies(
			viability_range = (4, 10),
            activation_range = (6, 10),
            maximum_growth_rate = 400u"m/Myr",
            extinction_coefficient = 0.1u"m^-1",
            saturation_intensity = 60u"W/m^2",
            diffusion_coefficient = 5000u"m"),

    	ALCAP2.Facies(
			viability_range = (4, 10),
            activation_range = (6, 10),
            maximum_growth_rate = 100u"m/Myr",
            extinction_coefficient = 0.005u"m^-1",
            saturation_intensity = 60u"W/m^2",
            diffusion_coefficient = 10000u"m")]

	input = ALCAP2.Input(
		tag = "ALCAP model",
		box = Common.Box{Shelf}(grid_size=(100, 50), phys_scale=150.0u"m"),
		time = TimeProperties(
			Δt = 200.0u"yr",
			steps = 5000,
			write_interval = 1),
		sea_level = t -> 4.0u"m" * sin(2π * t / 0.2u"Myr"),
		bedrock_elevation = (x, y) -> - x / 300.0,
		subsidence_rate = 50.0u"m/Myr",
		insolation = 400.0u"W/m^2",

    	disintegration_rate = 500.0u"m/Myr",
    	sediment_buffer_size = 50,
    	depositional_resolution = 0.5u"m",
		ca_interval = 1,

		facies = facies
	)

	H5Writer.run(Model{ALCAP2}, input, "../data/output/alcap1.h5")
end

# ╔═╡ 11fdbaa7-b4e1-4e41-9da2-c334f86e78ea
let
	header, data = read_slice("../data/output/alcap1.h5", :, 25)
	fig = Figure(size=(1000, 700))
	ax = Axis(fig[1,1])
	sediment_profile!(ax, header, data)
	fig
end

# ╔═╡ 4b8fa62e-ca50-4908-888e-6718ea6f5bf0
length(4:10)

# ╔═╡ 538cd660-aab5-409d-a508-89cfd62c128d
let
	header, data = read_slice("../data/output/alcap1.h5", :, 25)
	wheeler_diagram(header, data; smooth_size=(3,11), range=(-100.0u"m/Myr", 100.0u"m/Myr"))
end

# ╔═╡ Cell order:
# ╠═aa4723a2-8038-11ef-0a4b-a59fbd03c3ef
# ╠═19aca91a-f9af-4539-a91d-15923127cf9e
# ╠═20390248-fe1b-4dc2-bc28-4406fd7a91f3
# ╠═0c903628-0d63-4fcd-8e48-7351387f998b
# ╠═08c22253-a895-4d5f-93b1-74f628bd6b1b
# ╠═eb89211f-64d2-4cf7-a6ef-671dfabd4cc0
# ╠═bc3cc455-f7cb-4c57-9e7a-f78ae55afe14
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
# ╠═ae356c01-2520-4bdf-8555-8b465040e494
# ╠═c954e223-fb91-4f86-8199-f91f39156773
# ╟─e93b7ac8-dcf9-4ca1-8d1b-01d6dd139930
# ╠═bac48b9d-46d0-44a4-9e9f-7fa32b8fc75d
# ╠═79b3c99b-0691-405c-8ed8-1da6d25cb9a8
# ╠═67f405a6-1252-4fc8-9c70-0d97dd4437ae
# ╠═c8fbbf2a-de9e-4e88-b3bd-201091e8d869
# ╟─20901dd1-ce07-4dd0-aa1f-3673a8012da3
# ╠═a864d6ea-99bf-4513-9f44-f15b1d19d6af
# ╠═b523780a-7841-4f89-b144-f7b277bec831
# ╟─72e1d116-06f0-4623-b068-9dd32c6f0eba
# ╠═c3b7220f-49f4-42a2-bd2a-1c554e86d3c8
# ╠═11fdbaa7-b4e1-4e41-9da2-c334f86e78ea
# ╠═4b8fa62e-ca50-4908-888e-6718ea6f5bf0
# ╠═538cd660-aab5-409d-a508-89cfd62c128d
