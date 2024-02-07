### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ be3739c1-bbf3-4e48-8ced-f249197c2c98
using Pkg; Pkg.activate("../workenv")

# ╔═╡ ab41cc36-45d2-48eb-b549-edd7f49ece01
using Revise

# ╔═╡ 7f499c8a-e89f-499a-b19c-3ccd97a87ff8
using CarboKitten.Transport: deposit, Box, Particle

# ╔═╡ f168a0f8-0ef3-46aa-bac8-c82e3ebf6116
using CarboKitten.BoundaryTrait: Periodic

# ╔═╡ a0263acb-8f42-4b72-821f-c279d2042122
using GLMakie

# ╔═╡ 9e993f30-cf67-497e-b091-7edbe18bad61
const TestParticle = Particle{Nothing}

# ╔═╡ 14dac041-96e7-4746-b998-344ded3b84da
box = Box((16, 16), (x=1.0, y=1.0), 1.0/16.0)

# ╔═╡ 1614e9fe-a56f-4094-94e1-8eb296cc8024
function axes(box::Box)
	x_axis = 0:box.phys_scale:box.phys_size.x
	y_axis = 0:box.phys_scale:box.phys_size.y
	(x_axis, y_axis)
end

# ╔═╡ b011b065-3a36-4b50-9fce-a3faf6bb1f76
axes(box)

# ╔═╡ 1dabd7fb-4562-49c6-a8fa-980499d0737b
Periodic{2}

# ╔═╡ 71d4ab7b-c43c-4427-9ccc-3009c1257f60
particles = rand(Float64, 2, 42) |> eachcol .|> 
		(p -> TestParticle((x=p[1], y=p[2]), 1.0, 1.0, 1, nothing));

# ╔═╡ 43130df6-4c32-404a-9675-64ac1d058fad
density = let
	target = zeros(Float64, 1, 16, 16);
	particles .|> deposit(Periodic{2}, box, target)
	target
end;

# ╔═╡ fc2171b8-0ce1-4b99-8bbe-b4097814eb35
let
	fig = Figure()
	ax = Axis(fig[1, 1])
	heatmap!(ax, axes(box)..., density[1,:,:])
	scatter!(ax, [(p.position.x, p.position.y) for p in particles]; color=:red)
	fig
end

# ╔═╡ 0e4dd1a6-6b7e-439e-82ee-55d4e156d846
sum(density)

# ╔═╡ 720a1a2e-20d9-41d7-b680-5f94f876fa69
mod(1.0, 1.0)

# ╔═╡ Cell order:
# ╠═be3739c1-bbf3-4e48-8ced-f249197c2c98
# ╠═ab41cc36-45d2-48eb-b549-edd7f49ece01
# ╠═7f499c8a-e89f-499a-b19c-3ccd97a87ff8
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
