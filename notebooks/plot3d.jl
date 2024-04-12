### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 94b47a0a-f849-11ee-1d11-234cda7e097b
using Pkg; Pkg.activate("../workenv")

# ╔═╡ 519bba17-7a32-4309-a447-874d6e83a6ef
using Revise

# ╔═╡ d5511180-f859-4b83-9f48-f7e2d6b07242
using HDF5

# ╔═╡ 0c26747f-e423-488a-a648-f6eb1363bbaf
using GLMakie

# ╔═╡ ccbdb01c-c4ba-437e-8605-428059f8ade4
na = [CartesianIndex()]

# ╔═╡ 7d1b870a-739c-40c9-a0cc-ae44a04714a5
h = h5open("../data/caps-osc.h5", "r") do fid
	attr = HDF5.attributes(fid["input"])
	Δt = attr["delta_t"][]
	subsidence_rate = attr["subsidence_rate"][]
	t_end = fid["input/t"][end-1]
	total_subsidence = subsidence_rate * t_end
	total_sediment = sum(fid["sediment"][]; dims=3)
	initial_height = fid["input/height"][]
	center = div(size(total_sediment)[1], 2)
	n_t = size(total_sediment)[4]
	print(size(total_sediment))
	ev = cumsum(total_sediment[:,:,1,1:1000]; dims=3) .* Δt .- initial_height[na,:] .- total_subsidence
	# e1 = sum(total_sediment[:,:,:,1:(n_t÷4)]; dims=4) .* Δt .- initial_height[na,:] .- total_subsidence
	# e2 = sum(total_sediment[:,:,:,1:(n_t÷2)]; dims=4) .* Δt .- initial_height[na,:] .- total_subsidence
	# e4 = sum(total_sediment[:,:,:,:]; dims=4) .* Δt .- initial_height[na,:] .- total_subsidence
	t = fid["input/t"][]
	n_facies = size(fid["sediment"])[3]

	return ev
end

# ╔═╡ d9d83aa4-5c99-4295-a6d1-d3bc580fa443
let
	fig = Figure(size=(1000,800))
	ax = Axis3(fig[1,1]; aspect=(1, 2, 1), azimuth=0.1π)
	args = (:colormap => :Spectral, :colorrange => (-100, 0))
	x = 1:50 |> collect
	y = 1:100 |> collect
	surface!(ax, h[:,:,250]; args...)
	wireframe!(ax, x, y, h[x,y,250]; color=:black, alpha=0.1)
	surface!(ax, h[:,:,500]; args...)
	wireframe!(ax, x, y, h[x,y,500]; color=:black, alpha=0.1)
	surface!(ax, h[:,:,750]; args...)
	wireframe!(ax, x, y, h[x,y,750]; color=:black, alpha=0.1)
	surface!(ax, h[:,:,1000]; args...)
	wireframe!(ax, x, y, h[x,y,1000]; color=:black, alpha=0.1)
	save("../data/plot3d.png", fig; px_per_unit=2)
	fig
end

# ╔═╡ Cell order:
# ╠═94b47a0a-f849-11ee-1d11-234cda7e097b
# ╠═519bba17-7a32-4309-a447-874d6e83a6ef
# ╠═d5511180-f859-4b83-9f48-f7e2d6b07242
# ╠═ccbdb01c-c4ba-437e-8605-428059f8ade4
# ╠═7d1b870a-739c-40c9-a0cc-ae44a04714a5
# ╠═0c26747f-e423-488a-a648-f6eb1363bbaf
# ╠═d9d83aa4-5c99-4295-a6d1-d3bc580fa443
