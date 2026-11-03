# ~/~ begin <<docs/src/visualization/topography.md#ext/GlamourView.jl>>[init]
module GlamourView

import CarboKitten.Visualization: glamour_view!
using CarboKitten.Utility: in_units_of
using CarboKitten.Output.Storage: Header, DataVolume, Pack, PackVolume,
    sediment_thickness, x_axis, y_axis, time_axis, bedrock
using Makie
using Unitful

function glamour_view!(ax::Makie.Axis3, p::PackVolume; colormap=Reverse(:speed))
    x = x_axis(p) |> in_units_of(u"km")
    y = y_axis(p) |> in_units_of(u"km")
    xy_aspect = x[end] / y[end]

    ax.aspect = (xy_aspect, 1, 1)
    ax.azimuth = -π/3

    n_steps = length(time_axis(p))
    grid_size = (length(x), length(y))

    steps_between = 2
    selected_steps::Vector{Int} = if n_steps > (2 + steps_between)
        [1, ((1:steps_between) .* n_steps .÷ (steps_between + 1))..., n_steps]
    else
        collect(1:n_steps)
    end

    base = bedrock(p)
    result = Array{Float64, 3}(undef, grid_size..., length(selected_steps))
    st = sediment_thickness(p)
    for (i, j) in enumerate(selected_steps)
        result[:, :, i] = (st[:,:,j] .+ base) |> in_units_of(u"m")
    end

    surface!(ax, x, y, result[:,:,1];
    color=ones(grid_size),
    colormap=:grays)

    for s in eachslice(result[:,:,2:end-1], dims=3)
        surface!(ax, x, y, s;
        colormap=(colormap, 0.7))
    end

    surface!(ax, x, y, result[:,:,end]; colormap=colormap)
    lines!(ax, x, zeros(grid_size[1]), result[:, 1, end]; color=(:white, 0.5), linewidth=1)
    lines!(ax, fill(x[end], grid_size[2]), y, result[end, :, end]; color=(:white, 0.5), linewidth=1)
end

function glamour_view(p::PackVolume; colormap=Reverse(:speed))
    fig = Figure()
    ax = Axis3(fig[1, 1])
    glamour_view!(ax, p, colormap=colormap)
    return fig, ax
end

# Methods accepting a `Header`, `DataVolume` pair, kept for backward compatibility.
glamour_view!(ax::Makie.Axis3, header::Header, data::DataVolume; kwargs...) =
    glamour_view!(ax, Pack(header, data); kwargs...)
glamour_view(header::Header, data::DataVolume; kwargs...) =
    glamour_view(Pack(header, data); kwargs...)

end
# ~/~ end
