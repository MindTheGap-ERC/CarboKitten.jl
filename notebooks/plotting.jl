### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 4a12c6f6-425f-4cc7-baf1-bfb87766bcae
using Pkg; Pkg.activate("../workenv")

# ╔═╡ 1b1d9efa-1b42-475d-ad94-7c357521b95b
using CairoMakie

# ╔═╡ 4683a8f4-3ca9-4996-81b7-d96af3f53dd2
using HDF5

# ╔═╡ c1819c71-2c1b-4eab-bf0a-00d3ffbb683c
using GeometryBasics

# ╔═╡ b141d394-1290-4e69-b4da-a8980cb020c1


# ╔═╡ 62c10025-e6d0-42b8-83f0-ddf434555b34
CairoMakie.activate!()

# ╔═╡ 87dad1f7-b8fb-4e56-95ba-394d511914d3
let
	f = Figure()
	ax = Axis(f[1, 1])
	x = [1, 2, 1.5, 2]
	y = [1, 2, 0.5, 1.5]
	pts = Point{2,Float64}.(x, y)
	faces = [TriangleFace(3, 2, 1), TriangleFace(2, 3, 4)]
	m = GeometryBasics.normal_mesh(pts, faces)
	mesh!(ax, m, color=[1, 2, 3, 4])
	f
end

# ╔═╡ f827abe4-5091-410e-9f4a-894597584bbb


# ╔═╡ 76e71ad0-863d-4d7a-b0c2-c714c9d249db
h5open("../data/catp.h5", "r") do fid
	total_sediment = sum(fid["sediment"][]; dims=1)
	total_sediment
end


# ╔═╡ 8d4fc57a-05ed-45f5-834f-f597a478402d
elevation = h5open("../data/catp.h5","r") do fid
	attr = HDF5.attributes(fid["input"])
	Δt = attr["delta_t"][]
	subsidence_rate = attr["subsidence_rate"][]
	t_end = fid["input/t"][end-1]
	total_subsidence = subsidence_rate * t_end
	total_sediment = sum(fid["sediment"][]; dims=1)
	initial_height = fid["input/height"][:,:]
	elevation = sum(total_sediment; dims=4)[1,:,:,:] .* Δt .- initial_height .- total_subsidence
	t = fid["input/t"][]
	return elevation[:,:,1]
end

# ╔═╡ 96273124-cf66-4be7-91e7-a6e99ab1c331
elevation

# ╔═╡ 05a3d577-9a08-4c4b-88a1-0992e56e8d29
surface(elevation)

# ╔═╡ 06ca199e-642c-4d4e-a465-85589e268fe0
x, t, h, p = h5open("../data/catp.h5","r") do fid
	y = 30
	attr = HDF5.attributes(fid["input"])
	Δt = attr["delta_t"][]
	subsidence_rate = attr["subsidence_rate"][]
	t_end = fid["input/t"][end-1]
	total_subsidence = subsidence_rate * t_end
	total_sediment = sum(fid["sediment"][]; dims=1)
	initial_height = fid["input/height"][:,y]
	elevation = cumsum(total_sediment; dims=4)[1,:,y,:] .* Δt .- initial_height .- total_subsidence
	t = fid["input/t"][]
	return fid["input/x"][], [t; Δt*attr["time_steps"][]], hcat(.- initial_height .- total_subsidence, elevation), fid["sediment"][:,:,y,:]
end

# ╔═╡ f8969a09-426c-4670-8756-12802711ed33
x

# ╔═╡ bd05564a-3cf9-4c8d-aec6-84d9c9ad75c4
let
	pts = vec(Point{2,Float64}.(x, h[:,2:end]))
	c = vec(argmax(p; dims=1)[1,:,:] .|> (c -> c[1]))
	rect = Rect2(0.0, 0.0, 1.0, 1.0)
	m_tmp = GeometryBasics.mesh(Tesselation(rect, (100, 1000)))
	m = GeometryBasics.normal_mesh(pts, faces(m_tmp))

	f = Figure()
	ax = Axis(f[1, 1], xlabel="location", ylabel="depth", limits=((-12,100000), nothing))
	mesh!(ax, m, color=c, alpha=0.7)
	for idx in [1,501,1001]
		lines!(ax, x, h[:, idx], color=:black)
		text!(ax, -2.0, h[1, idx]; text="$(t[idx]) Myr", align=(:right, :center))
	end
	for idx in [250,750]
		lines!(ax, x, h[:, idx], color=:black, linewidth=0.5)
	end
	save("crosssection.png", f)
	f
end

# ╔═╡ 96163e97-df27-466f-b607-db888e9b83fa
clip(x) = x < 0 ? 0.0 : x

# ╔═╡ 9a07a7e2-242b-4979-904d-91fff447fcc5
begin
	q = rand(Float64, 10, 10) .- 0.5
	q = clip.(q)
end

# ╔═╡ 8f055c57-e0a8-4f42-bd48-aa8df927b571


# ╔═╡ Cell order:
# ╠═4a12c6f6-425f-4cc7-baf1-bfb87766bcae
# ╠═b141d394-1290-4e69-b4da-a8980cb020c1
# ╠═1b1d9efa-1b42-475d-ad94-7c357521b95b
# ╠═4683a8f4-3ca9-4996-81b7-d96af3f53dd2
# ╠═c1819c71-2c1b-4eab-bf0a-00d3ffbb683c
# ╠═62c10025-e6d0-42b8-83f0-ddf434555b34
# ╠═87dad1f7-b8fb-4e56-95ba-394d511914d3
# ╠═f827abe4-5091-410e-9f4a-894597584bbb
# ╠═76e71ad0-863d-4d7a-b0c2-c714c9d249db
# ╠═8d4fc57a-05ed-45f5-834f-f597a478402d
# ╠═96273124-cf66-4be7-91e7-a6e99ab1c331
# ╠═05a3d577-9a08-4c4b-88a1-0992e56e8d29
# ╠═06ca199e-642c-4d4e-a465-85589e268fe0
# ╠═f8969a09-426c-4670-8756-12802711ed33
# ╠═bd05564a-3cf9-4c8d-aec6-84d9c9ad75c4
# ╠═96163e97-df27-466f-b607-db888e9b83fa
# ╠═9a07a7e2-242b-4979-904d-91fff447fcc5
# ╠═8f055c57-e0a8-4f42-bd48-aa8df927b571
