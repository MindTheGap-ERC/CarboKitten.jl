# ~/~ begin <<docs/src/visualization.md#ext/WheelerDiagram.jl>>[init]
module WheelerDiagram

import CarboKitten.Visualization: wheeler_diagram, wheeler_diagram!
using CarboKitten.Export: Header, Data, DataSlice, read_data, read_slice
using CarboKitten.Utility: in_units_of
using Makie
using Unitful
using CarboKitten.BoundaryTrait
using CarboKitten.Stencil: convolution


const na = [CartesianIndex()]

elevation(h::Header, d::DataSlice) =
    let bl = h.bedrock_elevation[d.slice..., na],
        sr = h.axes.t[end] * h.subsidence_rate

        bl .+ d.sediment_elevation .- sr
    end

water_depth(header::Header, data::DataSlice) =
    let h = elevation(header, data),
        s = header.subsidence_rate .* (header.axes.t .- header.axes.t[end]),
        l = header.sea_level

        h .- (s.+l)[na, :]
    end

const Rate = typeof(1.0u"m/Myr")

function sediment_accumulation!(ax::Axis, header::Header, data::DataSlice;
    smooth_size::NTuple{2,Int}=(3, 11),
    colormap=Reverse(:curl),
    range::NTuple{2,Rate}=(-100.0u"m/Myr", 100.0u"m/Myr"))
    magnitude = sum(data.deposition .- data.disintegration; dims=1)[1, :, :] ./ (header.Δt * header.write_interval)
    blur = convolution(Shelf, ones(Float64, smooth_size...) ./ *(smooth_size...))
    wd = zeros(Float64, length(header.axes.x), length(header.axes.t))
    blur(water_depth(header, data) / u"m", wd)
    mag = zeros(Float64, length(header.axes.x), length(header.axes.t) - 1)
    blur(magnitude / u"m/Myr", mag)

    ax.ylabel = "time [Myr]"
    ax.xlabel = "position [km]"

    sa = heatmap!(ax, header.axes.x / u"km", header.axes.t / u"Myr", mag;
        colormap=colormap, colorrange=range ./ u"m/Myr")
    contour!(ax, header.axes.x / u"km", header.axes.t / u"Myr", wd;
        levels=[0], color=:red, linewidth=2, linestyle=:dash)
    return sa
end

function dominant_facies!(ax::Axis, header::Header, data::DataSlice;
    smooth_size::NTuple{2,Int}=(3, 11),
    colors=Makie.wong_colors())
    n_facies = size(data.production)[1]
    colormax(d) = getindex.(argmax(d; dims=1)[1, :, :], 1)

    dominant_facies = colormax(data.deposition)
    blur = convolution(Shelf, ones(Float64, smooth_size...) ./ *(smooth_size...))
    wd = zeros(Float64, length(header.axes.x), length(header.axes.t))
    blur(water_depth(header, data) / u"m", wd)

    ax.ylabel = "time [Myr]"
    ax.xlabel = "position [km]"

    ft = heatmap!(ax, header.axes.x / u"km", header.axes.t / u"Myr", dominant_facies;
        colormap=cgrad(colors[1:n_facies], n_facies, categorical=true),
        colorrange=(0.5, n_facies + 0.5))
    contourf!(ax, header.axes.x / u"km", header.axes.t / u"Myr", wd;
        levels=[0.0, 100.0], colormap=Reverse(:grays))
    contour!(ax, header.axes.x / u"km", header.axes.t / u"Myr", wd;
        levels=[0], color=:black, linewidth=2)
    return ft
end

function wheeler_diagram!(ax1::Axis, ax2::Axis, header::Header, data::DataSlice;
    smooth_size::NTuple{2,Int}=(3, 11),
    range::NTuple{2,Rate}=(-100.0u"m/Myr", 100.0u"m/Myr"))

    linkyaxes!(ax1, ax2)
    sa = sediment_accumulation!(ax1, header, data; smooth_size=smooth_size, range=range)
    ft = dominant_facies!(ax2, header, data; smooth_size=smooth_size)
    ax2.ylabel = ""

    return sa, ft
end

function wheeler_diagram(header::Header, data::DataSlice;
    smooth_size::NTuple{2,Int}=(3, 11),
    range::NTuple{2,Rate}=(-100.0u"m/Myr", 100.0u"m/Myr"))

    fig = Figure(size=(1000, 600))
    ax1 = Axis(fig[2, 1])
    ax2 = Axis(fig[2, 2])

    sa, ft = wheeler_diagram!(ax1, ax2, header, data; smooth_size=smooth_size, range=range)

    Colorbar(fig[1, 1], sa; vertical=false, label="sediment accumulation [m/Myr]")
    Colorbar(fig[1, 2], ft; vertical=false, ticks=1:3, label="dominant facies")
    fig
end

end
# ~/~ end