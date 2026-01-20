using CairoMakie
using CarboKitten.Visualization: sediment_profile, summary_plot, stratigraphic_column!
using CarboKitten.Export: read_slice, Header, DataSlice, read_column, DataColumn
using Unitful

function plot_sediment_profile(HDF5_file::String)
    header, data = read_slice(HDF5_file, :profile)
	fig = sediment_profile(header, data, show_unconformities = false)
    save("docs/src/_fig/dissolution_sediment_profile.png", fig)
end

const HDF5_file = "data/output/dissolution.h5"

header, data = read_slice(HDF5_file, :profile)

const location = 25

function plot_sediment_accumulation(header::Header, data::DataSlice, location::Int; ax::Axis)
    time_interval = (header.axes.t[end] - header.axes.t[1]) /
                     (size(data.sediment_thickness, 2) - 1)
    times_range = header.axes.t[1]:time_interval:header.axes.t[end]
    times = collect(times_range)
    times = Float64.(vec(times) ./ u"Myr")

    thickness = Float64.(vec(data.sediment_thickness[location, :]) ./ u"m")

    lines!(ax, times, thickness, color = :black)
end

function plot_stratigraphic_column(header::Header, data::DataSlice, location::Int; ax::Axis)
    fig = Figure()
    ax = Axis(fig[1, 1]; title="Stratigraphic Column at Location $(location)", ylabel="Depth [m]")
    
    column = data[location]
    
    stratigraphic_column!(ax, header, column; color = Makie.wong_colors())
    
    return fig
end


function plot_barrel(header::Header, data::DataSlice, location::Int)
    fig = Figure()

    ax1 = Axis(fig[1, 1]; title = "Sediment Accumulation at Location $(location)", xlabel = "Time (Myr)", ylabel = "Sediment Thickness (m)")
    plot_sediment_accumulation(header, data, location; ax = ax1)

    ax2 = Axis(fig[1, 2]; title = "Stratigraphic Column)", ylabel = "Depth (m)")
    colsize!(fig.layout, 2, Relative(0.15))
    
    linkaxes!(ax1, ax2)
    
    column = data[location]
    stratigraphic_column!(ax2, header, column; color = Makie.wong_colors())
    save("docs/src/_fig/dissolution_barrel_location_$(location).png", fig)
    return fig
end

plot_barrel(header, data, location)
plot_sediment_profile(HDF5_file)