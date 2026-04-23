

using Makie
import CarboKitten.Visualization:
     fence_plot, facies_colormap

export fence_plot

function fence_plot(
    grid::AbstractArray{<:Integer,3};
    topk = nothing,
    category_names::Vector{String},
    dx::Real,
    dz::Real,
    colormap = nothing,
    facies_colors = nothing,   # ← add this
    colorbar_label::String = "",
    ys = nothing,
    xs = nothing,
    VE::Real = 10,
    fence_alpha = 1.0,
    shading = false,
    frame_lw::Real = 1.2,
    frame_color = :black
)

if topk !== nothing
    nz = maximum(topk)
    grid = grid[:, :, 1:nz]
end

grid_plot = Float32.(grid)
grid_plot[grid_plot .== 0] .= NaN32

    Nx, Ny, Nz = size(grid)

    ys === nothing && (ys = round.(Int, range(0.1Ny, 0.9Ny; length=4)))
    xs === nothing && (xs = round.(Int, range(0.1Nx, 0.9Nx; length=4)))

    x_m = ((1:Nx) .- 0.5) .* dx
    y_m = ((1:Ny) .- 0.5) .* dx
    z   = ((1:Nz) .- 0.5) .* dz
    z_plot = z .* VE

    z1, z2 = z_plot[1], z_plot[end]
    ncat = length(category_names)


if colormap === nothing
    colormap = facies_colormap(
        ncat;
        facies_colors = facies_colors,
        include_nodeposit = false)
end
    fig = Figure(size=(1300, 800))
    ax = Axis3(fig[1,1],
        xlabel="x",
        ylabel="y",
        zlabel="z",
        aspect=:data
    )
    ax.zreversed = true

    Xxz = repeat(x_m, 1, Nz)
    Zxz = repeat(z_plot', Nx, 1)
    Yyz = repeat(y_m, 1, Nz)
    Zyz = repeat(z_plot', Ny, 1)

		function frame_y!(x1, x2, y0)
    kw = (; color=frame_color, linewidth=frame_lw, depth_shift=-1f-4)

    lines!(ax, [x1, x1], [y0, y0], [z1, z2]; kw...)
    lines!(ax, [x2, x2], [y0, y0], [z1, z2]; kw...)
    lines!(ax, [x1, x2], [y0, y0], [z2, z2]; kw...)
end


function frame_x!(x0, y1, y2)
    kw = (; color=frame_color, linewidth=frame_lw, depth_shift=-1f-4)

    lines!(ax, [x0, x0], [y1, y1], [z1, z2]; kw...)
    lines!(ax, [x0, x0], [y2, y2], [z1, z2]; kw...)
    lines!(ax, [x0, x0], [y1, y2], [z2, z2]; kw...)
end

    for y0 in ys
        y0m = y_m[y0]
        surface!(ax, Xxz, fill(y0m, Nx, Nz), Zxz;
            color = grid_plot[:, y0, :],
            colormap = colormap,
            colorrange = (0.5, ncat + 0.5),
            shading = false,
            alpha = fence_alpha
        )
        frame_y!(x_m[1], x_m[end], y0m)
    end

    for x0 in xs
        x0m = x_m[x0]
        surface!(ax, fill(x0m, Ny, Nz), Yyz, Zyz;
            color = grid_plot[x0, :, :],
            colormap = colormap,
            colorrange = (0.5, ncat + 0.5),
            shading = false,
            alpha = fence_alpha
        )
        frame_x!(x0m, y_m[1], y_m[end])
    end
	
	# Draw intersection lines between x and y panels
	function intersection_lines!(xs, ys)
		kw = (; color = frame_color, linewidth = frame_lw * 1.2, overdraw = true)

		for x0 in xs
		x0m = x_m[x0]
		for y0 in ys
			y0m = y_m[y0]

		lines!(ax,
		[x0m, x0m],
		[y0m, y0m],
		[z1, z2],
		color=:black,
		linewidth=2,
		depth_shift=-1f-4
	)
		end
	end
	end

	intersection_lines!(xs, ys)

    autolimits!(ax)
    lims = ax.finallimits[]
    zmin = lims.origin[3]
    zmax = lims.origin[3] + lims.widths[3]
    zt = collect(range(zmin, zmax; length=6))
    ax.zticks = (zt, string.(round.(zt ./ VE; digits=1)))

    Colorbar(fig[1,2],
    colormap = colormap,
    limits   = (0.5, ncat + 0.5),
    ticks    = (1:ncat, category_names),
    label    = colorbar_label
)

    return fig
end


