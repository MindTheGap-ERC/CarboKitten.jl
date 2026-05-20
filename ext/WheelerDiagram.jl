# ~/~ begin <<docs/src/visualization.md#ext/WheelerDiagram.jl>>[init]
module WheelerDiagram

import CarboKitten.Visualization: wheeler_diagram, wheeler_diagram!
using CarboKitten.Export: Header, Data, DataSlice, read_data, read_slice
using CarboKitten.Utility: in_units_of
using CarboKitten.Output.Abstract: stratigraphic_column
using Makie
using Unitful
using CarboKitten.BoundaryTrait
using CarboKitten.Stencil: convolution


const na = [CartesianIndex()]

elevation(h::Header, d::DataSlice) =
    let bl = h.initial_topography[d.slice..., na],
        sr = (h.axes.t[end] - h.axes.t[1]) * h.subsidence_rate

        bl .+ d.sediment_thickness .- sr
    end

water_depth(header::Header, data::DataSlice) =
    let h = elevation(header, data),
        wi = data.write_interval,
        s = header.subsidence_rate .* (header.axes.t[1:wi:end] .- header.axes.t[end]),
        l = header.sea_level[1:wi:end]

        h .- (s.+l)[na, :]
    end

const Rate = typeof(1.0u"m/Myr")

function sediment_accumulation!(ax::Axis, header::Header, data::DataSlice;
    smooth_size::NTuple{2,Int}=(3, 11),
    colormap=Reverse(:curl),
    range::NTuple{2,Rate}=(-100.0u"m/Myr", 100.0u"m/Myr"))
    wi = data.write_interval
    magnitude = sum(data.deposition .- data.disintegration; dims=1)[1, :, :] ./ (header.Δt * wi)
    blur = convolution(Shelf, ones(Float64, smooth_size...) ./ *(smooth_size...))
    wd = zeros(Float64, length(header.axes.x), length(header.axes.t[1:wi:end]))
    blur(water_depth(header, data) / u"m", wd)
    mag = zeros(Float64, length(header.axes.x), length(header.axes.t[1:wi:end]))
    blur(magnitude / u"m/Myr", mag)

    ax.ylabel = "time [Myr]"
    ax.xlabel = "position [km]"
    xkm = header.axes.x |> in_units_of(u"km")
    tmyr = header.axes.t[1:wi:end] |> in_units_of(u"Myr")

    sa = heatmap!(ax, xkm, tmyr, mag;
        colormap=colormap, colorrange=range ./ u"m/Myr")
    #contour!(ax, xkm, tmyr, wd;
    #    levels=[0], color=:red, linewidth=2, linestyle=:dash)
    return sa
end

function dominant_facies!(ax::Axis, header::Header, data::DataSlice;
    show::Symbol = :both,
    smooth_size::NTuple{2,Int} = (3, 11),
    colors = Makie.wong_colors())

    if show ∉ [:model, :preserved, :both]
        error("expected argument `show` to be one of `:model`, `:preserved`, `:both`; got $(show)")
    end

    prec = 10^-8 # prec. in m - below acc. is considered 0

    n_facies = size(data.production)[1]
    colormax(d) = getindex.(argmax(d; dims=1)[1, :, :], 1)

    wi = data.write_interval
    blur = convolution(Shelf, ones(Float64, smooth_size...) ./ *(smooth_size...))
    wd = zeros(Float64, length(header.axes.x), length(header.axes.t[1:wi:end]))
    blur(water_depth(header, data) / u"m", wd)

    ax.ylabel = "time [Myr]"
    ax.xlabel = "position [km]"

    xkm = header.axes.x |> in_units_of(u"km")
    tmyr = header.axes.t[1:wi:end] |> in_units_of(u"Myr")


    ft = if show == :model
        dominant_facies = colormax(data.deposition)
        dominant_facies = Matrix{Union{Missing, Int}}(dominant_facies)
        dominant_facies[ wd .> 0] .= missing
        heatmap!(ax, xkm, tmyr, dominant_facies;
            colormap=cgrad(colors[1:n_facies], n_facies, categorical=true),
            colorrange=(0.5, n_facies + 0.5),
            nan_color=:white)
    elseif show == :preserved
        sc = stratigraphic_column(data)
        dominant_facies = colormax(sc)
        dominant_facies = Matrix{Union{Missing, Int}}(dominant_facies)
        combined_acc = dropdims(sum(sc, dims = 1), dims = 1) |> in_units_of(u"m")
        dominant_facies[ combined_acc .< prec] .= missing
        heatmap!(ax, xkm, tmyr, dominant_facies;
            colormap=cgrad(colors[1:n_facies], n_facies, categorical=true),
            colorrange=(0.5, n_facies + 0.5),
            nan_color=:white)
    else
        sc = stratigraphic_column(data)
        dominant_facies_model = colormax(data.deposition)
        dominant_facies_model = Matrix{Union{Missing, Int}}(dominant_facies_model)
        dominant_facies_model[ wd .> 0] .= missing
        heatmap!(ax, xkm, tmyr, dominant_facies_model;
            colormap=cgrad(colors[1:n_facies], n_facies, categorical=true),
            colorrange=(0.5, n_facies + 0.5), alpha=0.3,
            nan_color=:white)
        dominant_facies_preserved = colormax(sc)
        dominant_facies_preserved = Matrix{Union{Missing, Int}}(dominant_facies_preserved)
        combined_acc = dropdims(sum(sc, dims = 1), dims = 1) |> in_units_of(u"m")
        dominant_facies_preserved[ combined_acc .< prec] .= missing
        heatmap!(ax, xkm, tmyr, dominant_facies_preserved;
            colormap=cgrad(colors[1:n_facies], n_facies, categorical=true),
            colorrange=(0.5, n_facies + 0.5),
            nan_color=:transparent)
    end

    #contourf!(ax, xkm, tmyr, wd;
    #    levels=[0.0, 10000.0], colormap=Reverse(:grays))
    #contour!(ax, xkm, tmyr, wd;
    #    levels=[0], color=:black, linewidth=2)
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
