### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ 54ca3310-b3a4-11ee-3914-fd0be29d888e
using Pkg; Pkg.activate("..")

# ╔═╡ 8039c056-4afb-4163-8157-37b8b611e91b
using GLMakie

# ╔═╡ 01779804-dbc3-47b1-bd66-28848b861d4a
using HDF5

# ╔═╡ 632e978b-1c52-444d-8248-0347123db467
using GeometryBasics

# ╔═╡ eab95149-3246-4223-a627-c7ca854b2f15
x = collect(-5:0.1:5)

# ╔═╡ ab28f836-90e0-407c-84ec-ec8cbf2ca2b4
y = cos.(x)

# ╔═╡ a7a4b5a5-cd12-4a1d-8de0-460bcdc4e8ff
let
	fig = Figure()
	ax = Axis(fig[1,1])
	lines!(ax, x, sin.(x))
	ax = Axis(fig[1,2])
	lines!(ax, x, cos.(x))
	fig
end

# ╔═╡ a4b9875a-588d-450d-b666-d9f6263732ff
# ╠═╡ disabled = true
#=╠═╡
sim_data = let
	h5open("../data/caps-osc.h5", "r") do fid
		attr = HDF5.attributes(fid["input"])
		Δt = attr["delta_t"][]
		initial_height = fid["input/height"][]
        subsidence_rate = attr["subsidence_rate"][]
        t_end = fid["input/t"][end-1]
        total_subsidence = subsidence_rate * t_end
		production_rate = sum(fid["sediment"][25,10,:,:]; dims=3)
		sediment = production_rate .* Δt
		age_depth = cumsum(sum(sediment; dims=1); dims=2) .- initial_height[10] .- total_subsidence

		(production_rate=production_rate, age_depth=age_depth)
	end
end
  ╠═╡ =#

# ╔═╡ bc0099e5-44a4-4bf4-a7db-46f47f4029f8
#=╠═╡
lines(collect(1:1000), sim_data.age_depth[1,:])
  ╠═╡ =#

# ╔═╡ 68448eae-b8d8-44ec-a39c-c06b244d1f38
function plot_crosssection(pos, datafile)
    # x: 1-d array with x-coordinates
    # t: 1-d array with time-coordinates (n_steps + 1)
    # h[x, t]: height fn, monotonic increasing in time
    # p[x, facies, t]: production rate
    # taken at y = y_max / 2, h[x, 1] is initial height
    n_facies, x, t, h, p = h5open(datafile,"r") do fid
        attr = HDF5.attributes(fid["input"])
        Δt = attr["delta_t"][]
        subsidence_rate = attr["subsidence_rate"][]
        t_end = fid["input/t"][end-1]
        total_subsidence = subsidence_rate * t_end
        total_sediment = sum(fid["sediment"][]; dims=3)
        initial_height = fid["input/height"][]
        center = div(size(total_sediment)[1], 2)
        elevation = cumsum(total_sediment; dims=4)[center,:,1,:] .* Δt .- initial_height .- total_subsidence
        t = fid["input/t"][]
        n_facies = size(fid["sediment"])[3]

        return n_facies,
               fid["input/x"][],
               [t; Δt*attr["time_steps"][]],
               hcat(.- initial_height .- total_subsidence, elevation),
               fid["sediment"][center,:,:,:]
    end

	pts = vec(Point{2,Float64}.(x, h[:,2:end]))
	c = vec(argmax(p; dims=2)[:,1,:] .|> (c -> c[2]))
	rect = Rect2(0.0, 0.0, 1.0, 1.0)
	m_tmp = GeometryBasics.mesh(Tesselation(rect, (100, 1000)))
	m = GeometryBasics.Mesh(pts, faces(m_tmp))

	# pts = vec(Point{2,Float64}.(x, h))
	# c = argmax(p; dims=2)[:,1,:] .|> (c -> c[2])
    # w = size(x)[1]

    # face(idx) = let k = idx[1] + idx[2]*w
    #     TriangleFace(k, k+1, k+1+w), TriangleFace(k+1+w, k+w, k)
    # end

	ax = Axis(pos, xlabel="location", ylabel="depth", limits=((-12,x[end]), nothing))
    # for f in 1:n_facies
    #     locs = CartesianIndices((size(x)[1], size(t)[1] - 1))[c .== f]
    #     triangles = collect(Iterators.flatten(face.(locs)))
    #     m = GeometryBasics.Mesh(pts, triangles)
    #     mesh!(ax, m)
    # end

	mesh!(ax, m, color=c, alpha=0.7)
	for idx in [1,501,1001]
		lines!(ax, x, h[:, idx], color=:black)
		text!(ax, -2.0, h[1, idx]; text="$(t[idx]) Myr", align=(:right, :center))
	end
	for idx in [250,750]
		lines!(ax, x, h[:, idx], color=:black, linewidth=0.5)
	end
end


# ╔═╡ 7d5ed268-1cb5-43fb-93ff-16da28769617
let
	fig = Figure()
	plot_crosssection(fig[1,1], "../data/caps-osc.h5")
	fig
end

# ╔═╡ 9fbb5dc1-2ba8-4576-bc32-4e72c448ecc9
sim_data2 = let
	h5open("../data/caps-osc.h5", "r") do fid
		attr = HDF5.attributes(fid["input"])
		Δt = attr["delta_t"][]
		initial_height = fid["input/height"][]
        subsidence_rate = attr["subsidence_rate"][]
        t_end = fid["input/t"][end-1]
        total_subsidence = subsidence_rate * t_end
		production_rate = sum(fid["sediment"][25,10,:,:]; dims=3)
		sediment = production_rate .* Δt
        erosion = sum(fid["erosion"][25,10,:]; dims=3)
		age_depth = cumsum(sum(sediment; dims=1); dims=2) .- cumsum(sum(erosion;
		dims=1);dims=2) .- initial_height[10] .- total_subsidence

		(production_rate=production_rate, age_depth=age_depth, erosion = erosion)
	end
end

# ╔═╡ d17fcbc8-ee50-4079-a23f-6ffc15412244
lines(collect(1:1000), sim_data2.age_depth[1,:])

# ╔═╡ Cell order:
# ╠═54ca3310-b3a4-11ee-3914-fd0be29d888e
# ╠═8039c056-4afb-4163-8157-37b8b611e91b
# ╠═01779804-dbc3-47b1-bd66-28848b861d4a
# ╠═632e978b-1c52-444d-8248-0347123db467
# ╠═eab95149-3246-4223-a627-c7ca854b2f15
# ╠═ab28f836-90e0-407c-84ec-ec8cbf2ca2b4
# ╠═a7a4b5a5-cd12-4a1d-8de0-460bcdc4e8ff
# ╠═a4b9875a-588d-450d-b666-d9f6263732ff
# ╠═bc0099e5-44a4-4bf4-a7db-46f47f4029f8
# ╠═68448eae-b8d8-44ec-a39c-c06b244d1f38
# ╠═7d5ed268-1cb5-43fb-93ff-16da28769617
# ╠═9fbb5dc1-2ba8-4576-bc32-4e72c448ecc9
# ╠═d17fcbc8-ee50-4079-a23f-6ffc15412244
