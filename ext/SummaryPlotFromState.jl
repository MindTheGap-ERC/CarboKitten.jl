module SummaryPlotFromState

using Makie
using Unitful
using CarboKitten.Utility: in_units_of
using CarboKitten.Export: Header
using CarboKitten.Export: dataslice_from_state_exact
using CarboKitten.Visualization:
    sediment_profile!, wheeler_diagram!, glamour_view!, production_curve!
using CarboKitten.Components.Production: production_rate

import CarboKitten.Visualization: summary_plot_from_state, facies_colormap



function summary_plot_from_state(
    header,
    state,
    input;
    wheeler_smooth = (1,1),
    show_unconformities = true,
    facies_colors = nothing,
    facies_names = nothing
)


    slice = dataslice_from_state_exact(state; y_index = 1, write_interval = 1)

    n_facies = length(input.facies)

    if facies_names === nothing
    facies_names = if all(hasproperty.(input.facies, :name))
        [getproperty(f, :name) for f in input.facies]
    else
        string.(1:length(input.facies))
    end
end

    fig = Figure(size = (1200, 1000), backgroundcolor = :gray80)

    # -------- Sediment profile --------
    ax1 = Axis(fig[1:2, 1:2])
    sediment_profile!(
    ax1, header, slice;
    show_unconformities = show_unconformities,
    show_coeval_lines = true,
    show_sealevel = true,
    facies_colors = facies_colors   # keep this
)
    # -------- Wheeler --------
    ax2 = Axis(fig[4,1])
    ax3 = Axis(fig[4,2])

    sm, df = wheeler_diagram!(
    ax2, ax3, header, slice;
    smooth_size = wheeler_smooth,
    facies_colors = facies_colors
)

    Colorbar(fig[3,1], sm;
    vertical = false,
    label = "sediment accumulation [m/Myr]"
)

Colorbar(fig[3,2], df;
    vertical = false,
    label = "facies",
    ticks = 0:n_facies
)
    # -------- Sea level --------
    ax4 = Axis(fig[4,3], title = "sea level", xlabel = "sea level [m]")
    lines!(
        ax4,
        header.sea_level |> in_units_of(u"m"),
        header.axes.t |> in_units_of(u"Myr")
    )

    # -------- 3D topography (compaction-aware) --------
    ax5 = Axis3(fig[1,3], title = "topography")
    glamour_view!(ax5, header, input, state)

    # -------- Production curves --------
ax6 = Axis(fig[2,3], title = "production curves")

production_curve!(
    ax6,
    input;
    facies_colors = facies_colors,
    max_depth = 500.0u"m",
    facies_names = facies_names
)

return fig   # ← add this line
end

end # module