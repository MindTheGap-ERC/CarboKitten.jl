# ~/~ begin <<docs/src/visualization.md#ext/GlamourView.jl>>[init]
module GlamourView

import CarboKitten.Visualization: glamour_view!
using CarboKitten.Utility: in_units_of
using CarboKitten.Export: read_header
using Makie
using HDF5
using Unitful

function glamour_view!(ax::Makie.Axis3, fid::HDF5.File; colormap=Reverse(:speed))
	header = read_header(fid)
	x = header.axes.x |> in_units_of(u"km")
	y = header.axes.y |> in_units_of(u"km")
	xy_aspect = x[end] / y[end]

	ax.aspect = (xy_aspect, 1, 1)
	ax.azimuth = -π/3

	n_steps = length(header.axes.t)
	grid_size = (length(x), length(y))
	steps_between = 2
	selected_steps = [1, ((1:steps_between) .* n_steps .÷ (steps_between + 1))..., n_steps]
	bedrock = header.bedrock_elevation .- header.axes.t[end] * header.subsidence_rate

	result = Array{Float64, 3}(undef, grid_size..., length(selected_steps))
	for (i, j) in enumerate(selected_steps)
		result[:, :, i] = fid["sediment_height"][:,:,j] .+ bedrock / u"m"
	end

	surface!(ax, x, y, result[:,:,1];
		color=ones(grid_size),
		colormap=:grays)

	for s in eachslice(result[:,:,2:end-1], dims=3)
		surface!(ax, x, y, s;
			colormap=(colormap, 0.7))		
	end

	surface!(ax, x, y, result[:,:,end];
		colormap=colormap)
	lines!(ax, x, zeros(grid_size[1]), result[:, 1, end]; color=(:white, 0.5), linewidth=1)
	lines!(ax, fill(x[end], grid_size[2]), y, result[end, :, end]; color=(:white, 0.5), linewidth=1)
end

end
# ~/~ end