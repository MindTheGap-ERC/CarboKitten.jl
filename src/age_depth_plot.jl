using HDF5
using GLMakie
using CSV
using DataFrames

function plot_age_depth(datafile,position::Int64)
	h5open(datafile, "r") do fid
		attr = HDF5.attributes(fid["input"])
		Δt = attr["delta_t"][]
		initial_height = fid["input/height"][]
        subsidence_rate = attr["subsidence_rate"][]
        t_end = fid["input/t"][end-1]
        total_subsidence = subsidence_rate * t_end
		production_rate = sum(fid["sediment"][25,position,:,:]; dims=1)
		sediment = production_rate .* Δt
        erosion = fid["erosion"][25,position,:] .* Δt
        print(size(erosion), "   ", size(sediment))
		age_depth = cumsum(sediment .- erosion';dims=2) .- initial_height[position] .- total_subsidence
		print(size(age_depth))
		(sediment=sediment, age_depth=age_depth,erosion=erosion,ssd=t_end)
	end
end



result = plot_age_depth("data/caps-test.h5",10)
result2 = plot_age_depth("data/caps-test.h5",80)
z = Figure()
ax = Axis(z[1,1],title = "Age_depth model (no denudation)", xlabel = "Age (kyr)", ylabel = "Depth (m)")#
z,l1 = plot!(ax,collect(1:1000),result.age_depth[:],color = :blue)
l2 = plot!(ax,collect(1:1000),result2.age_depth[:],color = :red)
z
save("data/withdenudation_miller.png",z)
df = DataFrame(Depth_shallow = result.age_depth[:],Depth_deep = result2.age_depth[:])
CSV.write("data/withdenudation.csv",df)
