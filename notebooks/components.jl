### A Pluto.jl notebook ###
# v0.19.45

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

# ╔═╡ 89854ba4-10b9-4250-a26a-676f4012a540
using CarboKitten.Components

# ╔═╡ 741cb70c-5335-453f-95a9-d7de05ddc297
using GLMakie

# ╔═╡ c742188b-caca-4f81-85ec-9b5c8712377b
@compose module BS92
  @mixin UniformProduction

  using ..Common
  using CSV
  using DataFrames
  using Interpolations
  using ..UniformProduction: uniform_production
  using ..TimeIntegration
  using ..WaterDepth

  function State(input::Input)
    ti_state = TimeIntegration.State(input)
    sediment_height = zeros(Height, input.box.grid_size...)
    return State(0, sediment_height)
  end

  function step(input::Input)
    τ = uniform_production(input)
    function (state::State)
	  prod = τ(state) .* input.time.Δt
      Δη = sum(prod; dims=1)[1,:,:]
      state.sediment_height .+= Δη
      state.step += 1
	  return prod
    end
  end

  function sealevel_curve()
       data = DataFrame(CSV.File("../data/bs92-sealevel-curve.csv"))
       linear_interpolation(data.time, data.depth)
  end

  struct Frame
	  deposition::Array{Amount, 3}
	  sediment_height::Array{Amount, 2}
  end

  function run(input::Input)
      step! = step(input)
      getwd = WaterDepth.water_depth(input)
      state = State(input)

      n_writes = input.time.steps ÷ input.time.write_interval
	  Channel{Frame}() do ch
	      for i = 1:n_writes
			  prod = zeros(Amount, n_facies(input), input.box.grid_size...)
	          for _ = 1:input.time.write_interval
				  prod .+= step!(state)
	          end
			  put!(ch, Frame(prod, copy(state.sediment_height)))
	      end
	  end
  end
end

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

# ╔═╡ 80508f8b-15a5-4b72-bbdc-49205e3102bd
let
	input = BS92.Input(
      box = Common.Box{Shelf}(grid_size=(100, 1), phys_scale=600.0u"m"),
      time = TimeProperties(
        Δt = 10.0u"yr",
        steps = 8000,
        write_interval = 100),
      sea_level = let sc = BS92.sealevel_curve()
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

    result = BS92.run(input)
    fig = Figure()
    ax = Axis(fig[1,1], xlabel="x (km)", ylabel="z (m)")
    x, y = box_axes(input.box)
    η0 = input.bedrock_elevation.(x, y')

    for l in result
        η = η0 .+ l.sediment_height
        lines!(ax, x |> in_units_of(u"km"), vec(η) |> in_units_of(u"m"), color=:steelblue4)
	end

  	fig
end

# ╔═╡ e93b7ac8-dcf9-4ca1-8d1b-01d6dd139930
md"""
# Multiple facies
"""

# ╔═╡ 79b3c99b-0691-405c-8ed8-1da6d25cb9a8


# ╔═╡ 20901dd1-ce07-4dd0-aa1f-3673a8012da3
md"""
# Production + CA
"""

# ╔═╡ 618b202b-dc3a-443e-aabb-5326a5983f79


# ╔═╡ Cell order:
# ╠═aa4723a2-8038-11ef-0a4b-a59fbd03c3ef
# ╠═19aca91a-f9af-4539-a91d-15923127cf9e
# ╠═20390248-fe1b-4dc2-bc28-4406fd7a91f3
# ╠═0c903628-0d63-4fcd-8e48-7351387f998b
# ╠═08c22253-a895-4d5f-93b1-74f628bd6b1b
# ╠═eb89211f-64d2-4cf7-a6ef-671dfabd4cc0
# ╠═89854ba4-10b9-4250-a26a-676f4012a540
# ╠═741cb70c-5335-453f-95a9-d7de05ddc297
# ╟─3adfb906-aef1-4262-89e0-ba6c318c1ebc
# ╠═e486467a-af3c-4735-b336-d989fc31f355
# ╠═54b3f849-2b3d-4c8b-b932-a4980e54c354
# ╟─f6025f38-9417-43bd-ae48-f0528eb7025c
# ╠═c742188b-caca-4f81-85ec-9b5c8712377b
# ╠═80508f8b-15a5-4b72-bbdc-49205e3102bd
# ╟─e93b7ac8-dcf9-4ca1-8d1b-01d6dd139930
# ╠═79b3c99b-0691-405c-8ed8-1da6d25cb9a8
# ╟─20901dd1-ce07-4dd0-aa1f-3673a8012da3
# ╠═618b202b-dc3a-443e-aabb-5326a5983f79
