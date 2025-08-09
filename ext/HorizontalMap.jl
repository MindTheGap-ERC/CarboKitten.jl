# ~/~ begin <<docs/src/visualization.md#ext/HorizontalMap.jl>>[init]
module HorizontalMap

using Makie
using Unitful

using CarboKitten.Visualization
using CarboKitten.Export

"""
    plot_dominant_facies_heatmap(data_array, time_step; title_prefix="Dominant Facies")

Create a heatmap showing the dominant facies at a specific time step.

# Arguments
- `data_array`: 4D array with dimensions (facies, x, y, time), usually data.deposition but works for others
- `time_step`: Integer index for the time step to visualize
- `title_prefix`: Prefix for the plot title

# Returns
- Makie Figure object containing the heatmap
"""

function plot_dominant_facies_heatmap(data_array, time_step; 
                                     title_prefix="Dominant Facies")
    time_slice = data_array[:, :, :, time_step]
    
    dominant_facies_indices = argmax(time_slice; dims=1)[1, :, :]
    dominant_facies = getindex.(dominant_facies_indices, 1)
    
    n_facies = size(data_array, 1)
    
    fig = Figure()
    ax = Axis(fig[1, 1], 
              title="$title_prefix at Time Step $time_step",
              xlabel="Distance [grid cell units]", 
              ylabel="Distance [grid cell units]")
    
    hm = heatmap!(ax, dominant_facies, 
                  colormap=cgrad(Makie.wong_colors()[1:n_facies], n_facies, categorical=true),
                  colorrange=(0.5, n_facies + 0.5))
    
    Colorbar(fig[1, 2], hm, 
             label="Facies Index", 
             ticks=1:n_facies)
    
    return fig
end


end
# ~/~ end
