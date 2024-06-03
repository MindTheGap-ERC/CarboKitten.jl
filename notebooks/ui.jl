### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 9fbb6b40-f1e4-11ee-36ea-85e4b86183c8
using Pkg; Pkg.activate("../workenv")

# ╔═╡ 6060e8f9-4a91-4559-a194-335b8d4a7d54
using CarboKitten.BoundaryTrait: Shelf

# ╔═╡ 48a5634a-480d-433e-9bf6-614bdfca75b0
using CarboKitten.Config: Box, TimeProperties

# ╔═╡ 3647417e-b0b5-4dce-b43c-63d9bfb5a12f
using CarboKitten.Burgess2013: production_rate, Facies

# ╔═╡ f74b79b1-b715-42b5-bcf8-560950b77676
using CarboKitten.CaProd

# ╔═╡ 908447dc-fb5b-4aa6-95ff-11ca69eec349
using GLMakie

# ╔═╡ 2976428b-06f9-4447-9d29-04b0486bec33
using Unitful

# ╔═╡ 884dc010-4c7a-41dc-9939-aaee75d060f0
using Observables

# ╔═╡ 4db60e9d-62ed-49c5-9bbb-2586cec2d079
begin
	const PERIOD = 200.0u"kyr"
	const AMPLITUDE = 4.0u"m"
	const FACIES = [
	    Facies(viability_range = (4, 10),
	           activation_range = (6, 10),
	           maximum_growth_rate = 500u"m/Myr",
	           extinction_coefficient = 0.8u"m^-1",
	           saturation_intensity = 60u"W/m^2"),
	
	    Facies(viability_range = (4, 10),
	           activation_range = (6, 10),
	           maximum_growth_rate = 400u"m/Myr",
	           extinction_coefficient = 0.1u"m^-1",
	           saturation_intensity = 60u"W/m^2"),
	
	    Facies(viability_range = (4, 10),
	           activation_range = (6, 10),
	           maximum_growth_rate = 100u"m/Myr",
	           extinction_coefficient = 0.005u"m^-1",
	           saturation_intensity = 60u"W/m^2")
	]

	const INPUT = CaProd.Input(
	  box = Box{Shelf}(
	    grid_size = (100, 50),
	    phys_scale = 0.15u"km"
	  ),
	  time = TimeProperties(
	    Δt = 0.001u"Myr",
	    steps = 1000,
	    write_interval = 1
	  ),
	  sea_level = t -> AMPLITUDE * sin(2π * t / PERIOD), 
	  subsidence_rate=50.0u"m/Myr",
	  initial_depth=x -> x / 300.0,
	  facies=FACIES,
	  insolation=400.0u"W/m^2"
	)
end

# ╔═╡ 3d69675d-03b1-420c-9f23-83cf81aca304
15_000/50

# ╔═╡ 4391904d-29cf-4202-9698-39104c265dd5
function plot_facies_production(input; loc = nothing)
	fig, loc = isnothing(loc) ? let fig = Figure(); (fig, fig[1, 1]) end : (nothing, loc)
	ax = Axis(loc, title="production at $(sprint(show, INPUT.insolation; context=:fancy_exponent=>true))", xlabel="production (m/Myr)", ylabel="depth (m)", yreversed=true)
	for f in input.facies
		depth = (0.1:0.1:50.0)u"m"
		prod = [production_rate(input.insolation, f, d) for d in depth]
		lines!(ax, prod / u"m/Myr", depth / u"m")
		
	end
	fig
end

# ╔═╡ e26b52ea-ac60-4486-8171-b5b072850f3c
plot_facies_production(INPUT)

# ╔═╡ c5a5a050-4dcb-4cfc-8451-6f78d79f8fb0
height = Observable{Matrix{Float64}}(zeros(Float64, INPUT.box.grid_size...))

# ╔═╡ 4ca647bb-0b47-4e4e-8871-d8d48c9b20ca
function run_model(input)
    Channel{CaProd.Frame}() do ch
        s = CaProd.initial_state(input)
        p = CaProd.propagator(input)
        u = CaProd.updater(input)

        while true
            Δ = p(s)
            put!(ch, Δ)
            u(s, Δ)
			# height[] = s.height / u"m"
        end
    end
end

# ╔═╡ ff34f5c8-88e0-4f54-af66-cda661a77d32
for f in Iterators.take(run_model(INPUT), 1000)
end

# ╔═╡ 09b87b3a-ad32-49c3-b69f-14d2c8e147ee
"production at $(sprint(show, INPUT.insolation; context=:fancy_exponent=>true))"

# ╔═╡ Cell order:
# ╠═9fbb6b40-f1e4-11ee-36ea-85e4b86183c8
# ╠═6060e8f9-4a91-4559-a194-335b8d4a7d54
# ╠═48a5634a-480d-433e-9bf6-614bdfca75b0
# ╠═3647417e-b0b5-4dce-b43c-63d9bfb5a12f
# ╠═f74b79b1-b715-42b5-bcf8-560950b77676
# ╠═908447dc-fb5b-4aa6-95ff-11ca69eec349
# ╠═2976428b-06f9-4447-9d29-04b0486bec33
# ╠═884dc010-4c7a-41dc-9939-aaee75d060f0
# ╠═4db60e9d-62ed-49c5-9bbb-2586cec2d079
# ╠═3d69675d-03b1-420c-9f23-83cf81aca304
# ╠═4391904d-29cf-4202-9698-39104c265dd5
# ╠═e26b52ea-ac60-4486-8171-b5b072850f3c
# ╠═c5a5a050-4dcb-4cfc-8451-6f78d79f8fb0
# ╠═4ca647bb-0b47-4e4e-8871-d8d48c9b20ca
# ╠═ff34f5c8-88e0-4f54-af66-cda661a77d32
# ╠═09b87b3a-ad32-49c3-b69f-14d2c8e147ee
