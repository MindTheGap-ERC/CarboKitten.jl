# ~/~ begin <<docs/src/visualization.md#ext/GlamourView.jl>>[init]
module GlamourView

import CarboKitten.Visualization: glamour_view!
using CarboKitten.Utility: in_units_of
using CarboKitten.Export: Header, DataVolume
using Makie
using HDF5
using Unitful

function glamour_view!(ax::Makie.Axis3, header::Header, data::DataVolume; colormap=Reverse(:speed))
    x = header.axes.x[data.slice[1]] |> in_units_of(u"km")
    y = header.axes.y[data.slice[2]] |> in_units_of(u"km")
	xy_aspect = x[end] / y[end]

	ax.aspect = (xy_aspect, 1, 1)
	ax.azimuth = -π/3

    n_steps = size(data.sediment_thickness)[3]
	grid_size = (length(x), length(y))
	steps_between = 2
	selected_steps = [1, ((1:steps_between) .* n_steps .÷ (steps_between + 1))..., n_steps]
	bedrock = header.initial_topography .- header.axes.t[end] * header.subsidence_rate

	result = Array{Float64, 3}(undef, grid_size..., length(selected_steps))
	for (i, j) in enumerate(selected_steps)
        result[:, :, i] = (data.sediment_thickness[:,:,j] .+ bedrock) |> in_units_of(u"m")
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

function glamour_view(header::Header, data::DataVolume; colormap=Reverse(:speed))
    fig = Figure()
    ax = Axis3(fig[1, 1])
    glamour_view!(ax, header, data, colormap=colormap)
    return fig, ax
end

end
# ~/~ end
