### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ a3dbe990-0cb1-11f1-a216-b9296323da7e
using Pkg; Pkg.activate("../workenv")

# ╔═╡ 1771b3bc-8058-443b-a769-f9bdad27ea01
using CarboKitten

# ╔═╡ ff22def0-dd36-4ae4-b982-b46aca9bd23f
using CarboKitten.Models: ALCAP as M

# ╔═╡ 05810aae-abd3-40c7-acae-c6fea42d2398
using GLMakie

# ╔═╡ 78816f97-7742-4571-a2b1-0c2fbabe3842
include("transport/diffusivity_estimation.jl")

# ╔═╡ 39dc0261-4c4b-4513-b0cb-4fc879336ea6
box = CarboKitten.Box{Periodic{2}}(grid_size=(100, 1), phys_scale=50.0u"m")

# ╔═╡ d9a525fe-29c1-42a8-841a-26d88d395796
begin
	t_end = 1.0u"Myr"
	Δt = 50.0u"yr"
	t_steps = t_end / Δt |> ceil |> Int
	peak_centre = box.phys_scale * box.grid_size[1] ÷ 2
	peak_width = 200.0u"m"
	peak_height = 10.0u"m"
	write_interval = max(1, t_steps ÷ 1000)
end

# ╔═╡ 81ab5596-8f07-4573-8088-bc4940d192f5
facies1 = M.Facies(
	diffusion_coefficient = 10.0u"m/yr",
	initial_sediment = (x, _) -> peak_height * exp(-(x - peak_centre)^2/(2 * peak_width^2)),
	active = false,
)

# ╔═╡ f4a4c3ba-8e17-459b-8e85-039c98a50056
input = M.Input(
	facies = [facies1],
	box = box,
	time = TimeProperties(Δt=Δt, steps=t_steps),
	output = Dict(:profile => OutputSpec(slice=(:, 1), write_interval=write_interval)),
	
	cementation_time = 100.0u"yr",
	disintegration_rate = 50.0u"m/Myr",
	subsidence_rate = 0.0u"m/Myr",
	sea_level = _ -> 0.0u"m",
	initial_topography = (_, _) -> -100.0u"m",
	insolation = 0.0u"W/m^2",
	
	transport_solver = Val{:forward_euler},
	depositional_resolution = 1.0u"km",
	sediment_buffer_size = 2,
)

# ╔═╡ f6e2fac9-bed0-4ab6-9e2d-3b24804dcc45
result = run_model(Model{M}, input, MemoryOutput(input))

# ╔═╡ 9811cb99-80f3-46ff-bc85-3b9b928f8b88
result.data_slices[:profile].sediment_thickness[:,1]

# ╔═╡ d2e0c3b8-5cd0-4dca-845e-e7a5b238f440
let fig = Figure()
	ax = Axis(fig[1, 1])
	x = result.header.axes.x |> in_units_of(u"km")
	y = result.data_slices[:profile].sediment_thickness |> in_units_of(u"m")
	for i in [1, 100, 1000]
		lines!(ax, x, y[:,i])
	end
	fig
end

# ╔═╡ fcd1c3ad-b3f9-4111-9887-94f01e996075
DiffusivityEstimation.estimate_diffusivity(result)

# ╔═╡ Cell order:
# ╠═a3dbe990-0cb1-11f1-a216-b9296323da7e
# ╠═1771b3bc-8058-443b-a769-f9bdad27ea01
# ╠═ff22def0-dd36-4ae4-b982-b46aca9bd23f
# ╠═05810aae-abd3-40c7-acae-c6fea42d2398
# ╠═39dc0261-4c4b-4513-b0cb-4fc879336ea6
# ╠═d9a525fe-29c1-42a8-841a-26d88d395796
# ╠═81ab5596-8f07-4573-8088-bc4940d192f5
# ╠═f4a4c3ba-8e17-459b-8e85-039c98a50056
# ╠═f6e2fac9-bed0-4ab6-9e2d-3b24804dcc45
# ╠═9811cb99-80f3-46ff-bc85-3b9b928f8b88
# ╠═d2e0c3b8-5cd0-4dca-845e-e7a5b238f440
# ╠═78816f97-7742-4571-a2b1-0c2fbabe3842
# ╠═fcd1c3ad-b3f9-4111-9887-94f01e996075
