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


# ╔═╡ 87dad1f7-b8fb-4e56-95ba-394d511914d3
let
	f = Figure()
	ax = Axis(f[1, 1])
	x = [1, 2, 1.5, 2]
	y = [1, 2, 0.5, 1.5]
	pts = Point{2,Float64}.(x, y)
	faces = [TriangleFace(3, 2, 1), TriangleFace(3, 2, 4)]
	m = GeometryBasics.Mesh(pts, faces)
	mesh!(ax, m, color=[1, 2, 3, 4])
	f
end

# ╔═╡ 76e71ad0-863d-4d7a-b0c2-c714c9d249db
h5open("../data/catp.h5", "r") do fid
	total_sediment = sum(fid["trans_sediment"][]; dims=3)
	total_sediment
end


# ╔═╡ 06ca199e-642c-4d4e-a465-85589e268fe0
x, t, h, p = h5open("../data/caps-osc.h5","r") do fid
	attr = HDF5.attributes(fid["input"])
	Δt = attr["delta_t"][]
	subsidence_rate = attr["subsidence_rate"][]
	t_end = fid["input/t"][end-1]
	total_subsidence = subsidence_rate * t_end
	total_sediment = sum(fid["sediment"][]; dims=3)
	initial_height = fid["input/height"][]
	elevation = cumsum(total_sediment; dims=4)[25,:,1,:] .* Δt .- initial_height .- total_subsidence
	t = fid["input/t"][]
	return fid["input/x"][], [t; Δt*attr["time_steps"][]], hcat(.- initial_height .- total_subsidence, elevation), fid["sediment"][25,:,:,:]
end

# ╔═╡ f8969a09-426c-4670-8756-12802711ed33


# ╔═╡ bd05564a-3cf9-4c8d-aec6-84d9c9ad75c4
let
	pts = vec(Point{2,Float64}.(x, h[:,2:end]))
	c = vec(argmax(p; dims=2)[:,1,:] .|> (c -> c[2]))
	rect = Rect2(0.0, 0.0, 1.0, 1.0)
	m_tmp = GeometryBasics.mesh(Tesselation(rect, (100, 1000)))
	m = GeometryBasics.Mesh(pts, faces(m_tmp))

	f = Figure()
	ax = Axis(f[1, 1], xlabel="location", ylabel="depth", limits=((-12,100), nothing))
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

# ╔═╡ 62bb362a-bd7e-4a4d-b9cf-4ff0dd3ea9f7
t

# ╔═╡ Cell order:
# ╠═4a12c6f6-425f-4cc7-baf1-bfb87766bcae
# ╠═b141d394-1290-4e69-b4da-a8980cb020c1
# ╠═1b1d9efa-1b42-475d-ad94-7c357521b95b
# ╠═4683a8f4-3ca9-4996-81b7-d96af3f53dd2
# ╠═c1819c71-2c1b-4eab-bf0a-00d3ffbb683c
# ╠═87dad1f7-b8fb-4e56-95ba-394d511914d3
# ╠═76e71ad0-863d-4d7a-b0c2-c714c9d249db
# ╠═06ca199e-642c-4d4e-a465-85589e268fe0
# ╠═f8969a09-426c-4670-8756-12802711ed33
# ╠═bd05564a-3cf9-4c8d-aec6-84d9c9ad75c4
# ╠═62bb362a-bd7e-4a4d-b9cf-4ff0dd3ea9f7
