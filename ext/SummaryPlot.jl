# ~/~ begin <<docs/src/visualization.md#ext/SummaryPlot.jl>>[init]
module SummaryPlot

using CarboKitten.Visualization
import CarboKitten.Visualization: summary_plot
using CarboKitten.Export: read_header, read_slice
using CarboKitten.Utility: in_units_of
using HDF5
using Unitful
using Makie

summary_plot(filename::AbstractString; kwargs...) = h5open(fid->summary_plot(fid; kwargs...), filename, "r")

function summary_plot(fid::HDF5.File; wheeler_smooth=(1, 1))
	header = read_header(fid)
    y_slice = div(length(header.axes.y), 2) + 1
    max_depth = minimum(header.bedrock_elevation)
	data = read_slice(fid, :, y_slice)

	n_facies = size(data.production)[1]
	fig = Figure(size=(1200, 1000), backgroundcolor=:gray80)

	ax1 = Axis(fig[1:2,1:2])
	sediment_profile!(ax1, header, data)

	ax2 = Axis(fig[4,1])
	ax3 = Axis(fig[4,2])
	sm, df = wheeler_diagram!(ax2, ax3, header, data; smooth_size=wheeler_smooth)
	Colorbar(fig[3,1], sm; vertical=false, label="sediment accumulation [m/Myr]")
	Colorbar(fig[3,2], df; vertical=false, label="dominant facies", ticks=1:n_facies)

	ax4 = Axis(fig[4,3], title="sealevel curve", xlabel="sealevel [m]",
        limits=(nothing, (header.axes.t[1] |> in_units_of(u"Myr"),
						  header.axes.t[end] |> in_units_of(u"Myr"))))
	lines!(ax4, header.sea_level |> in_units_of(u"m"), header.axes.t |> in_units_of(u"Myr"))

	ax5 = Axis(fig[2,3])
	production_curve!(ax5, fid["input"], max_depth=max_depth)

	linkyaxes!(ax2, ax3, ax4)

	ax = Axis3(fig[1, 3]; zlabel="depth [m]", xlabel="x [km]", ylabel="y [km]")
	glamour_view!(ax, fid)

	fig
end

end
# ~/~ end