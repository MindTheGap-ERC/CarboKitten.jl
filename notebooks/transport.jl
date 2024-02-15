### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ be3739c1-bbf3-4e48-8ced-f249197c2c98
using Pkg; Pkg.activate("../workenv")

# ╔═╡ ab41cc36-45d2-48eb-b549-edd7f49ece01
using Revise

# ╔═╡ 7f499c8a-e89f-499a-b19c-3ccd97a87ff8
using CarboKitten.Transport: deposit, Box, Particle, interpolate

# ╔═╡ bd30f1b8-891f-476b-8635-3f41eee710e3
using CarboKitten.Vectors

# ╔═╡ f168a0f8-0ef3-46aa-bac8-c82e3ebf6116
using CarboKitten.BoundaryTrait: Periodic, Reflected

# ╔═╡ a0263acb-8f42-4b72-821f-c279d2042122
using GLMakie

# ╔═╡ 6dde775c-00e3-4ddb-b1d7-62e40ec3df39
using CarboKitten.CATP: submarine_transport, State, ProductFrame, Input, Facies, stress, particles

# ╔═╡ f802a537-563e-4587-9e13-c3cb557b2947
using CarboKitten.Transport

# ╔═╡ 9e993f30-cf67-497e-b091-7edbe18bad61
const TestParticle = Particle{Nothing}

# ╔═╡ 14dac041-96e7-4746-b998-344ded3b84da
box16 = Box((16, 16), (x=1.0, y=1.0), 1.0/16.0)

# ╔═╡ 1614e9fe-a56f-4094-94e1-8eb296cc8024
function axes(box::Box)
	x_axis = 0:box.phys_scale:box.phys_size.x
	y_axis = 0:box.phys_scale:box.phys_size.y
	(x_axis, y_axis)
end

# ╔═╡ 1dabd7fb-4562-49c6-a8fa-980499d0737b
Periodic{2}

# ╔═╡ 71d4ab7b-c43c-4427-9ccc-3009c1257f60
test_particles = rand(Float64, 2, 42) |> eachcol .|> 
		(p -> TestParticle((x=p[1], y=p[2]), 1.0, 1.0, 1, nothing));

# ╔═╡ 43130df6-4c32-404a-9675-64ac1d058fad
density = let
	target = zeros(Float64, 1, 16, 16);
	test_particles .|> deposit(Periodic{2}, box16, target)
	target
end;

# ╔═╡ fc2171b8-0ce1-4b99-8bbe-b4097814eb35
let
	fig = Figure()
	ax = Axis(fig[1, 1])
	heatmap!(ax, axes(box16)..., density[1,:,:])
	scatter!(ax, [(p.position.x, p.position.y) for p in test_particles]; color=:red)
	fig
end

# ╔═╡ 0e4dd1a6-6b7e-439e-82ee-55d4e156d846
sum(density)

# ╔═╡ 720a1a2e-20d9-41d7-b680-5f94f876fa69
mod(1.0, 1.0)

# ╔═╡ c1a68718-1cbf-4abc-ad15-ee533ea51cd1
h = Float64[0 0; 1 1]

# ╔═╡ 0debe4b0-dcfe-4ec7-bb36-265fbd30bd03
box = Box((2, 2), (x=1.0, y=1.0), 0.5)

# ╔═╡ b011b065-3a36-4b50-9fce-a3faf6bb1f76
axes(box)

# ╔═╡ a6ae3820-5808-4e34-a090-f1b14a08f3f7
f = interpolate(Periodic{2}, box, h)

# ╔═╡ 862c53ce-1600-4348-9174-36949a72e19a
f((x=0.49, y=0.5))

# ╔═╡ 88421113-33ee-4d36-98ab-0ff2e558d361
(9, 7) .* 2

# ╔═╡ 4b6950f8-e02f-4149-8892-666e8f380eae
md"""
## Testing transport

To test transport, we set a height map and specify a production.
"""

# ╔═╡ a2b9f2f3-d213-4a42-8ecd-48a0a29d70c9
test_profile(x, y) = -x + 2.0 * exp(-((y-2.5)^2+(x-2.5)^2)/1.00)

# ╔═╡ 94dd0ded-69a2-4d13-b5a4-074d496669e1
input = Input(
	sea_level = t -> 0.0,
	subsidence_rate = 0.0,
	initial_depth = test_profile,
	grid_size = (50, 50),
	boundary = Reflected{2},
	phys_scale = 5.0/50,
	Δt = 1.0,
	time_steps = 100,
	write_interval = 1,
	facies = [ Facies((4, 10), (6, 10), 500.0, 0.8, 300, 1.0, 1.0, 7.4) ],
	insolation = 2000.0,
	Δz = 1.0,
	buffer_depth = 50,
	disintegration_rate = z -> 0.1,
	wave_shear_stress = z -> (x=0.0, y=0.0),
	g = 9.8,
	transport_subsample = 1
)

# ╔═╡ b9d09778-ccc6-4124-a027-51c775dc5171
begin
	x = collect(1:input.grid_size[1]) .* input.phys_scale
	y = collect(1:input.grid_size[2]) .* input.phys_scale
end

# ╔═╡ 3da02eb7-5eea-468a-bb3a-253c751b001b
height = test_profile.(x, y')

# ╔═╡ 6b332977-79cd-4d4b-ad00-30c00f2fab39
surface(x, x, height)

# ╔═╡ 91fc7690-c9c7-4a93-86d9-7e38e642ac34
sediment = zeros(Float64, input.grid_size..., input.buffer_depth, 1);

# ╔═╡ 50e6c1ba-c380-48b3-a0b5-b08461bf915b
state = State(0.0, height, sediment);

# ╔═╡ 0e11dbe1-db11-40e9-9dd3-a8faf1161aa7
production = ProductFrame(ones(Float64, 1, input.grid_size...) .* 0.1)

# ╔═╡ 0ad81d9c-0db5-4994-a5c1-9c72817e96d0
function submarine_transport_2(input::Input)
  # This function does not modify the state, rather it transports the
  # sediments in a given product frame and gives a new product frame
  # as output.
  box = Transport.Box(input.grid_size, (x=input.grid_size[1] * input.phys_scale, y=input.grid_size[2] * input.phys_scale), input.phys_scale)
  function (s::State, Δ::ProductFrame)  # -> ProductFrame
    output = zeros(Float64, size(Δ.production)...)
    transport = Transport.transport(input.boundary, box, stress(input, s))
    deposit = Transport.deposit(input.boundary, box, output)

    for p in particles(input, Δ)
      p |> transport |> deposit
    end

    return ProductFrame(output)
  end
end

# ╔═╡ e2cffe29-153a-48c5-948a-536a6cbf4894
tr = submarine_transport_2(input)

# ╔═╡ da298bd9-2f89-4929-bc5d-0c669edcf5b8
td = tr(state, production)

# ╔═╡ ff79db72-5086-45dc-975d-2423a0757049
let
	fig = Figure()
	ax = Axis(fig[1,1])
	co = heatmap!(ax, td.production[1,:,:])
	Colorbar(fig[1,2], co)
	fig
end

# ╔═╡ 982c8248-a6e2-4ce4-a20c-88f884dd9ef9
[(p.position.x, p.position.y) for p in particles(input, production)]
# stress(input, state)

# ╔═╡ d8a15e90-10d8-47c0-916c-cd4ddf8af2a6
points, stresses = let σ = stress(input, state)
	ps = collect(particles(input, production))
	points = [Point2(p.position.x, p.position.y) for p in ps]
	stresses = [let s = σ(p); Vec2f(s.x, s.y) end for p in ps]
	points, stresses
end

# ╔═╡ 5434f84b-4bc2-4d71-bad0-d66d7ccbd75c
strength = reshape(abs.(stress(input, state).(particles(input, production))), input.grid_size)

# ╔═╡ 68063b97-21a3-4413-ba31-26e9786d2739
arrows(points, stresses; lengthscale=0.01, arrowsize=5)

# ╔═╡ 29ab45a7-37de-4c71-b2c7-5696ecace703
let
	fig = Figure()
	ax = Axis(fig[1, 1])
	co = contourf!(ax, x, y, strength; levels=10)
	Colorbar(fig[1,2], co)
	fig
end

# ╔═╡ 343b01bf-114b-4371-82eb-e0fe11a6e4ba
function make_box(input)
    phys_size = (x = input.grid_size[1] * input.phys_scale,
                 y = input.grid_size[2] * input.phys_scale)
    box = Box(input.grid_size, phys_size, input.phys_scale)
	box
end

# ╔═╡ 0d32ed73-52a8-4a3d-a0a4-b262f4ad2022
make_box(input)

# ╔═╡ ffda1357-bdad-4249-83fa-d0d9c5d5a6e6
sum(strength) / length(strength)

# ╔═╡ Cell order:
# ╠═be3739c1-bbf3-4e48-8ced-f249197c2c98
# ╠═ab41cc36-45d2-48eb-b549-edd7f49ece01
# ╠═7f499c8a-e89f-499a-b19c-3ccd97a87ff8
# ╠═bd30f1b8-891f-476b-8635-3f41eee710e3
# ╠═f168a0f8-0ef3-46aa-bac8-c82e3ebf6116
# ╠═9e993f30-cf67-497e-b091-7edbe18bad61
# ╠═14dac041-96e7-4746-b998-344ded3b84da
# ╠═1614e9fe-a56f-4094-94e1-8eb296cc8024
# ╠═b011b065-3a36-4b50-9fce-a3faf6bb1f76
# ╠═1dabd7fb-4562-49c6-a8fa-980499d0737b
# ╠═71d4ab7b-c43c-4427-9ccc-3009c1257f60
# ╠═43130df6-4c32-404a-9675-64ac1d058fad
# ╠═a0263acb-8f42-4b72-821f-c279d2042122
# ╠═fc2171b8-0ce1-4b99-8bbe-b4097814eb35
# ╠═0e4dd1a6-6b7e-439e-82ee-55d4e156d846
# ╠═720a1a2e-20d9-41d7-b680-5f94f876fa69
# ╠═c1a68718-1cbf-4abc-ad15-ee533ea51cd1
# ╠═0debe4b0-dcfe-4ec7-bb36-265fbd30bd03
# ╠═a6ae3820-5808-4e34-a090-f1b14a08f3f7
# ╠═862c53ce-1600-4348-9174-36949a72e19a
# ╠═88421113-33ee-4d36-98ab-0ff2e558d361
# ╟─4b6950f8-e02f-4149-8892-666e8f380eae
# ╠═6dde775c-00e3-4ddb-b1d7-62e40ec3df39
# ╠═f802a537-563e-4587-9e13-c3cb557b2947
# ╠═94dd0ded-69a2-4d13-b5a4-074d496669e1
# ╠═a2b9f2f3-d213-4a42-8ecd-48a0a29d70c9
# ╠═b9d09778-ccc6-4124-a027-51c775dc5171
# ╠═3da02eb7-5eea-468a-bb3a-253c751b001b
# ╠═6b332977-79cd-4d4b-ad00-30c00f2fab39
# ╠═91fc7690-c9c7-4a93-86d9-7e38e642ac34
# ╠═50e6c1ba-c380-48b3-a0b5-b08461bf915b
# ╠═0e11dbe1-db11-40e9-9dd3-a8faf1161aa7
# ╠═0ad81d9c-0db5-4994-a5c1-9c72817e96d0
# ╠═e2cffe29-153a-48c5-948a-536a6cbf4894
# ╠═da298bd9-2f89-4929-bc5d-0c669edcf5b8
# ╠═ff79db72-5086-45dc-975d-2423a0757049
# ╠═982c8248-a6e2-4ce4-a20c-88f884dd9ef9
# ╠═d8a15e90-10d8-47c0-916c-cd4ddf8af2a6
# ╠═5434f84b-4bc2-4d71-bad0-d66d7ccbd75c
# ╠═68063b97-21a3-4413-ba31-26e9786d2739
# ╠═29ab45a7-37de-4c71-b2c7-5696ecace703
# ╠═343b01bf-114b-4371-82eb-e0fe11a6e4ba
# ╠═0d32ed73-52a8-4a3d-a0a4-b262f4ad2022
# ╠═ffda1357-bdad-4249-83fa-d0d9c5d5a6e6
