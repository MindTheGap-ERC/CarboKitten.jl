using Makie
import CarboKitten.Visualization: map_view, facies_colormap

export map_view

function map_view(
    grid::AbstractArray{<:Integer,3};
    z0::Int,
    category_names::Vector{String},
    category_colors = nothing,
    colorbar_label::String = "",
)

    Nx, Ny, Nz = size(grid)
    @assert 1 ≤ z0 ≤ Nz "z0 must be in 1:$Nz"

    slice_xy = grid[:, :, z0]

slice_plot = Float32.(slice_xy)
slice_plot[slice_plot .== 0] .= NaN32

    ncat = length(category_names)

    # Default categorical colormap
    cmap = facies_colormap(
    ncat;
    facies_colors = category_colors,
    include_nodeposit = false
)

    fig = Figure(size = (900, 700))
    ax  = Axis(fig[1,1],
        xlabel = "x index",
        ylabel = "y index",
        title  = "Map view (z = $z0)"
    )

    heatmap!(ax,
    slice_plot';
        colormap   = cmap,
        colorrange = (0.5, ncat + 0.5)
    )

    Colorbar(fig[1,2],
        colormap = cmap,
        limits = (0.5, ncat + 0.5),
ticks  = (1:ncat, category_names),
        label    = colorbar_label
    )

    return fig
end
