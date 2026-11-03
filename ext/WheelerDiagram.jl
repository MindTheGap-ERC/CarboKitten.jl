# ~/~ begin <<docs/src/visualization.md#ext/WheelerDiagram.jl>>[init]
module WheelerDiagram

import CarboKitten.Visualization: wheeler_diagram, wheeler_diagram!, sediment_accumulation!, dominant_facies!
using CarboKitten.Utility: in_units_of
using CarboKitten.Output.Storage: Header, DataSlice, Pack, PackSlice,
    stratigraphic_column, water_depth, deposition, net_deposition, n_facies,
    time_axis, write_time_step, section_axis
using Makie
using Unitful
using CarboKitten.BoundaryTrait
using CarboKitten.Stencil: convolution

const Rate = typeof(1.0u"m/Myr")

function sediment_accumulation!(ax::Axis, p::PackSlice;
    smooth_size::NTuple{2,Int}=(3, 11),
    colormap=Reverse(:curl),
    range::NTuple{2,Rate}=(-100.0u"m/Myr", 100.0u"m/Myr"))
    magnitude = dropdims(sum(net_deposition(p); dims=1), dims=1) ./ write_time_step(p)
    blur = convolution(Shelf, ones(Float64, smooth_size...) ./ *(smooth_size...))
    mag = zeros(Float64, size(magnitude))
    blur(magnitude / u"m/Myr", mag)

    ax.ylabel = "time [Myr]"
    ax.xlabel = "position [km]"
    xkm = section_axis(p) |> in_units_of(u"km")
    tmyr = time_axis(p) |> in_units_of(u"Myr")

    sa = heatmap!(ax, xkm, tmyr, mag;
        colormap=colormap, colorrange=range ./ u"m/Myr")
    return sa
end

function dominant_facies!(ax::Axis, p::PackSlice;
    show::Symbol = :both,
    smooth_size::NTuple{2,Int} = (3, 11),
    colors = Makie.wong_colors())

    if show ∉ [:model, :preserved, :both]
        error("expected argument `show` to be one of `:model`, `:preserved`, `:both`; got $(show)")
    end

    prec = 10^-8 # prec. in m - below acc. is considered 0

    n_f = n_facies(p)
    colormax(d) = getindex.(argmax(d; dims=1)[1, :, :], 1)

    blur = convolution(Shelf, ones(Float64, smooth_size...) ./ *(smooth_size...))
    wd_raw = water_depth(p) / u"m"
    wd = zeros(Float64, size(wd_raw))
    blur(wd_raw, wd)

    ax.ylabel = "time [Myr]"
    ax.xlabel = "position [km]"

    xkm = section_axis(p) |> in_units_of(u"km")
    tmyr = time_axis(p) |> in_units_of(u"Myr")

    ft = if show == :model
        dominant_facies = colormax(deposition(p))
        dominant_facies = Matrix{Union{Missing, Int}}(dominant_facies)
        dominant_facies[ wd .> 0] .= missing
        heatmap!(ax, xkm, tmyr, dominant_facies;
            colormap=cgrad(colors[1:n_f], n_f, categorical=true),
            colorrange=(0.5, n_f + 0.5),
            nan_color=:white)
    elseif show == :preserved
        sc = stratigraphic_column(p)
        dominant_facies = colormax(sc)
        dominant_facies = Matrix{Union{Missing, Int}}(dominant_facies)
        combined_acc = dropdims(sum(sc, dims = 1), dims = 1) |> in_units_of(u"m")
        dominant_facies[ combined_acc .< prec] .= missing
        heatmap!(ax, xkm, tmyr, dominant_facies;
            colormap=cgrad(colors[1:n_f], n_f, categorical=true),
            colorrange=(0.5, n_f + 0.5),
            nan_color=:white)
    else
        sc = stratigraphic_column(p)
        dominant_facies_model = colormax(deposition(p))
        dominant_facies_model = Matrix{Union{Missing, Int}}(dominant_facies_model)
        dominant_facies_model[ wd .> 0] .= missing
        heatmap!(ax, xkm, tmyr, dominant_facies_model;
            colormap=cgrad(colors[1:n_f], n_f, categorical=true),
            colorrange=(0.5, n_f + 0.5), alpha=0.3,
            nan_color=:white)
        dominant_facies_preserved = colormax(sc)
        dominant_facies_preserved = Matrix{Union{Missing, Int}}(dominant_facies_preserved)
        combined_acc = dropdims(sum(sc, dims = 1), dims = 1) |> in_units_of(u"m")
        dominant_facies_preserved[ combined_acc .< prec] .= missing
        heatmap!(ax, xkm, tmyr, dominant_facies_preserved;
            colormap=cgrad(colors[1:n_f], n_f, categorical=true),
            colorrange=(0.5, n_f + 0.5),
            nan_color=:transparent)
    end

    return ft
end

function wheeler_diagram!(ax1::Axis, ax2::Axis, p::PackSlice;
    smooth_size::NTuple{2,Int}=(3, 11),
    range::NTuple{2,Rate}=(-100.0u"m/Myr", 100.0u"m/Myr"))

    linkyaxes!(ax1, ax2)
    sa = sediment_accumulation!(ax1, p; smooth_size=smooth_size, range=range)
    ft = dominant_facies!(ax2, p; smooth_size=smooth_size)
    ax2.ylabel = ""

    return sa, ft
end

function wheeler_diagram(p::PackSlice;
    smooth_size::NTuple{2,Int}=(3, 11),
    range::NTuple{2,Rate}=(-100.0u"m/Myr", 100.0u"m/Myr"))

    fig = Figure(size=(1000, 600))
    ax1 = Axis(fig[2, 1])
    ax2 = Axis(fig[2, 2])

    sa, ft = wheeler_diagram!(ax1, ax2, p; smooth_size=smooth_size, range=range)

    Colorbar(fig[1, 1], sa; vertical=false, label="sediment accumulation [m/Myr]")
    Colorbar(fig[1, 2], ft; vertical=false, ticks=1:n_facies(p), label="dominant facies")
    fig
end

# Methods accepting a `Header`, `DataSlice` pair, kept for backward compatibility.
sediment_accumulation!(ax::Axis, header::Header, data::DataSlice; kwargs...) =
    sediment_accumulation!(ax, Pack(header, data); kwargs...)
dominant_facies!(ax::Axis, header::Header, data::DataSlice; kwargs...) =
    dominant_facies!(ax, Pack(header, data); kwargs...)
wheeler_diagram!(ax1::Axis, ax2::Axis, header::Header, data::DataSlice; kwargs...) =
    wheeler_diagram!(ax1, ax2, Pack(header, data); kwargs...)
wheeler_diagram(header::Header, data::DataSlice; kwargs...) =
    wheeler_diagram(Pack(header, data); kwargs...)

end
# ~/~ end
