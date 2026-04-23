using Makie
import CarboKitten.Visualization: extract_column

export extract_column

function extract_column(
    grid::AbstractArray{<:Integer,3},
    x::Int,
    y::Int;
    topk = nothing,
    category_names::Vector{String},
    category_colors = nothing,
    dz::Real,
    title_str::String = "",
    colorbar_label::String = ""
)


if topk !== nothing
    nz = maximum(topk)
    grid = grid[:, :, 1:nz]
end
    Nx, Ny, Nz = size(grid)

    @assert 1 ≤ x ≤ Nx
    @assert 1 ≤ y ≤ Ny

    ncat = length(category_names)

    # Default categorical colormap if none provided
    if category_colors === nothing
        cmap = Makie.cgrad(:Set3, ncat; categorical=true)
    else
        cmap = Makie.cgrad(Makie.to_color.(category_colors), ncat; categorical=true)
    end

    # Extract SINGLE voxel column
    column_codes = Float32.(grid[x, y, 1:Nz])

    zmid = ((1:Nz) .- 0.5) .* dz
    xgrid = [0.0, 1.0]
    img = repeat(reshape(column_codes, 1, :), 2, 1)

    fig = Figure(size=(360, 850))
    ax = Axis(fig[1,1];
        ylabel = "Depth",
        title  = isempty(title_str) ?
            "Borehole (x=$x, y=$y)" :
            title_str,
        yreversed = true,
        xticksvisible = false,
        xticklabelsvisible = false
    )

    heatmap!(ax, xgrid, zmid, img;
        colormap = cmap,
        colorrange = (-0.5f0, ncat - 0.5f0)
    )

    Colorbar(fig[1,2];
        colormap = cmap,
        limits = (-0.5, ncat - 0.5),
        ticks = (0:ncat-1, category_names),
        label = colorbar_label
    )

    return fig
end
