### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 3d000bc6-d6e8-11ee-039b-01742fb58212
using Pkg; Pkg.activate("../workenv")

# ╔═╡ 9c3a775e-046f-42b5-bf01-3f60c6dff053
using Revise

# ╔═╡ 628068ae-5acc-486c-a2ef-a1021fe1bdd7
using CarboKitten.BoundaryTrait

# ╔═╡ 2ce65a61-8703-42d2-8fc1-9208060e110d
using CarboKitten.Transport: interpolate, Box

# ╔═╡ e568b8ab-064d-48de-8cec-62aa957fbafd


# ╔═╡ b09002aa-924e-4685-8bb6-3730e11ba7ef
begin
	function axes(box::Box)
	    x = collect((0:box.grid_size[1]-1) .* box.phys_scale)
		y = collect((0:box.grid_size[2]-1) .* box.phys_scale)'
		return x, y
	end
	
	make_vec2(x, y) = (x=x, y=y)
	phys_grid(box::Box) = make_vec2.(axes(box)...)
end

# ╔═╡ 41c385bf-fe2c-48ab-8ce8-7dce38a9f9b1
box = Box((3, 2), (x=3.0, y=2.0), 1.0)

# ╔═╡ a27ad6b6-060d-4623-8b2a-24cd2830cbe7
values = collect([2.0 1.0 0.0; 1.0 2.0 -1.0]')

# ╔═╡ f5475bce-4170-499d-85e3-e0f8867bddad
interpolate(Shelf, box, values, (x=2.5, y=0.5))

# ╔═╡ Cell order:
# ╠═e568b8ab-064d-48de-8cec-62aa957fbafd
# ╠═3d000bc6-d6e8-11ee-039b-01742fb58212
# ╠═9c3a775e-046f-42b5-bf01-3f60c6dff053
# ╠═628068ae-5acc-486c-a2ef-a1021fe1bdd7
# ╠═2ce65a61-8703-42d2-8fc1-9208060e110d
# ╠═b09002aa-924e-4685-8bb6-3730e11ba7ef
# ╠═41c385bf-fe2c-48ab-8ce8-7dce38a9f9b1
# ╠═a27ad6b6-060d-4623-8b2a-24cd2830cbe7
# ╠═f5475bce-4170-499d-85e3-e0f8867bddad
