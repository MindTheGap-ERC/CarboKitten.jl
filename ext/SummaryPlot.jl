# ~/~ begin <<docs/src/visualization/overview.md#ext/SummaryPlot.jl>>[init]
module SummaryPlot

using CarboKitten.Visualization
import CarboKitten.Visualization: summary_plot
using CarboKitten.Output.Abstract: AbstractBundle, load, load_volume, load_slice
using CarboKitten.Output.Storage: header, data_volumes, data_slices, n_facies,
    time_axis, y_axis, sea_level
using CarboKitten.Output.H5Reader: H5Bundle, read_header
using CarboKitten.Output.MemoryWriter: MemoryOutput
using CarboKitten.Utility: in_units_of
using HDF5
using Unitful
using Makie

summary_plot(filename::AbstractString; kwargs...) = load(b -> summary_plot(b; kwargs...), filename)
summary_plot(fid::HDF5.File; kwargs...) = summary_plot(H5Bundle(fid, read_header(fid)); kwargs...)
summary_plot(output::MemoryOutput; kwargs...) = summary_plot(load(output); kwargs...)

# The production curve is read from the `input` group of an HDF5 file.
# Other bundles don't store this information.
_production_curve!(fig, bundle::H5Bundle) =
    production_curve!(Axis(fig[2,3]), bundle.fid["input"])
_production_curve!(fig, bundle::AbstractBundle) = nothing

function summary_plot(bundle::AbstractBundle; wheeler_smooth=(1, 1), show_unconformities=true)
    hdr = header(bundle)
    volumes = sort(data_volumes(bundle))
    slices = sort(data_slices(bundle))

    if isempty(slices) && isempty(volumes)
        @warn "No volume data or slice data stored. Cannot produce summary view."
        return nothing
    end

    fig = Figure(size=(1200, 1000), backgroundcolor=:gray80)

    volume_data = if isempty(volumes)
        @warn "No volume data stored, skipping topographic plots."
        nothing
    else
        if length(volumes) > 1
            @warn "Multiple volume data sets, picking first one."
        end

        volume_data = load_volume(bundle, volumes[1])
        ax = Axis3(fig[1, 3]; title="topography", zlabel="depth [m]", xlabel="x [km]", ylabel="y [km]")
        glamour_view!(ax, volume_data)
        volume_data
    end

    section_data = if isempty(slices)
        @warn "No profile data slice stored, taking section of volume data along x-axis."

        y_slice = div(length(y_axis(volume_data)), 2) + 1
        volume_data[:, y_slice]
    else
        load_slice(bundle, slices[1])
    end

    n_f = n_facies(section_data)

    ax1 = Axis(fig[1:2,1:2])
    sediment_profile!(ax1, section_data; show_unconformities = show_unconformities)
    axislegend(ax1; merge=true, backgroundcolor=:gray80)

    ax2 = Axis(fig[4,1])
    ax3 = Axis(fig[4,2])
    sm, df = wheeler_diagram!(ax2, ax3, section_data; smooth_size=wheeler_smooth)
    Colorbar(fig[3,1], sm; vertical=false, label="sedimentation rate [m/Myr]")
    Colorbar(fig[3,2], df; vertical=false, label="dominant facies", ticks=1:n_f)
    ax4 = Axis(fig[4,3], title="sealevel curve", xlabel="sealevel [m]",
               limits=(nothing, (first(time_axis(section_data)) |> in_units_of(u"Myr"),
                                 last(time_axis(section_data)) |> in_units_of(u"Myr"))))
    lines!(ax4, sea_level(hdr) |> in_units_of(u"m"), time_axis(hdr) |> in_units_of(u"Myr"))

    _production_curve!(fig, bundle)

    linkyaxes!(ax2, ax3, ax4)

    fig
end

end
# ~/~ end
