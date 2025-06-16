# ~/~ begin <<docs/src/visualization.md#ext/SummaryPlot.jl>>[init]
module SummaryPlot

using CarboKitten.Visualization
import CarboKitten.Visualization: summary_plot
using CarboKitten.Export: read_header, read_volume, read_slice, group_datasets
using CarboKitten.Utility: in_units_of
using HDF5
using Unitful
using Makie

summary_plot(filename::AbstractString; kwargs...) = h5open(fid->summary_plot(fid; kwargs...), filename, "r")

function summary_plot(fid::HDF5.File; wheeler_smooth=(1, 1))
    header = read_header(fid)
    data_groups = group_datasets(fid)

    if length(data_groups[:slice]) == 0 && length(data_groups[:volume]) == 0
        @warn "No volume data or slice data stored. Cannot produce summary view."
        return nothing
    end

    fig = Figure(size=(1200, 1000), backgroundcolor=:gray80)

    volume_data = if length(data_groups[:volume]) == 0
        @warn "No volume data stored, skipping topographic plots."
        nothing
    else
        if length(data_groups[:volume]) > 1
            @warn "Multiple volume data sets, picking first one."
        end

        volume_data = read_volume(fid[data_groups[:volume][1]])
        ax = Axis3(fig[1, 3]; zlabel="depth [m]", xlabel="x [km]", ylabel="y [km]")
        glamour_view!(ax, header, volume_data)
        volume_data
    end

    section_data = if length(data_groups[:slice]) == 0
        @warn "No profile data slice stored, taking section of volume data along x-axis."

        y_slice = div(size(volume_data.sediment_thickness)[2], 2) + 1
        volume_data[:, y_slice]
    else
        read_slice(fid[data_groups[:slice][1]])
    end

    n_facies = size(section_data.production)[1]

    ax1 = Axis(fig[1:2,1:2])
    sediment_profile!(ax1, header, section_data)

    ax2 = Axis(fig[4,1])
    ax3 = Axis(fig[4,2])
    sm, df = wheeler_diagram!(ax2, ax3, header, section_data; smooth_size=wheeler_smooth)
    Colorbar(fig[3,1], sm; vertical=false, label="sedimentation rate [m/Myr]")
    Colorbar(fig[3,2], df; vertical=false, label="dominant facies", ticks=1:n_facies)

    ax4 = Axis(fig[4,3], title="sealevel curve", xlabel="sealevel [m]",
               limits=(nothing, (header.axes.t[1] |> in_units_of(u"Myr"),
                                 header.axes.t[end] |> in_units_of(u"Myr"))))
    lines!(ax4, header.sea_level |> in_units_of(u"m"), header.axes.t |> in_units_of(u"Myr"))

    ax5 = Axis(fig[2,3])
    max_depth = minimum(header.initial_topography)
    production_curve!(ax5, fid["input"], max_depth=max_depth)

    linkyaxes!(ax2, ax3, ax4)

    fig
end

end
# ~/~ end
