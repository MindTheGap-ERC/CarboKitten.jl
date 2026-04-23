module WheelerDiagram

import CarboKitten.Visualization: wheeler_diagram, wheeler_diagram!, facies_colormap
using CarboKitten.Export: Header, DataSlice
using CarboKitten.Utility: in_units_of
using Makie
using Unitful
using CarboKitten.BoundaryTrait
using CarboKitten.Stencil: convolution

const Rate = typeof(1.0u"m/Myr")
const na = [CartesianIndex()]

# -----------------------------
# SAME SHAPES AS HDF5 EXPORT
# -----------------------------


water_depth(header::Header, data::DataSlice) =
    let base = header.initial_topography[data.slice..., na],
        subs = data.cumulative_subsidence,
        sea  = header.sea_level[1:size(subs,2)]

        base .+ data.sediment_thickness .- subs .- sea'
    end


function sediment_accumulation!(ax::Axis, header::Header, data::DataSlice;
    smooth_size=(3,11),
    colormap=Reverse(:curl),
    range=(-100.0u"m/Myr", 100.0u"m/Myr")
)
    wi = data.write_interval

    # Nt−1
    dep = data.deposition
    dis = data.disintegration
    mag_raw = sum(dep .- dis; dims=1)[1,:,:] ./ (header.Δt * wi)

    blur = convolution(Shelf, ones(Float64, smooth_size...) ./ *(smooth_size...))

    Nt1 = size(dep,3)
    raw_wd = water_depth(header, data) / u"m"
    raw_wd = raw_wd[:,1:Nt1]              # ensure Nt−1

    wd = zeros(Float64, size(raw_wd))
    blur(raw_wd, wd)

    mag = zeros(Float64, size(dep,2), Nt1)
    blur(mag_raw / u"m/Myr", mag)

    xkm  = header.axes.x |> in_units_of(u"km")
    tmyr = header.axes.t[1:Nt1] |> in_units_of(u"Myr")   # Nt−1 ticks

    heatmap!(ax, xkm, tmyr, mag; colormap, colorrange = range ./ u"m/Myr")
end


function dominant_facies!(ax::Axis, header::Header, data::DataSlice;
    smooth_size=(3,11),
    facies_colors=nothing,
    nondep_eps=0.0u"m"
)
    dep = data.deposition
    dom = getindex.(argmax(dep; dims=1)[1,:,:], 1)
    tot = sum(dep; dims=1)[1,:,:]

    domf = Float64.(dom)
    domf[tot .<= nondep_eps] .= NaN

    blur = convolution(Shelf, ones(Float64, smooth_size...) ./ *(smooth_size...))

    Nt1 = size(dep,3)
    out = zeros(Float64, size(dep,2), Nt1)
    blur(domf, out)

    out[out .< 0.5] .= NaN

    xkm  = header.axes.x |> in_units_of(u"km")
    tmyr = header.axes.t[1:Nt1] |> in_units_of(u"Myr")

    n_f = size(dep,1)
cmap = facies_colormap(n_f; facies_colors=facies_colors, include_nodeposit=true)

heatmap!(ax, xkm, tmyr, out;
    colormap = cmap,
    colorrange = (0.5f0, n_f + 0.5f0)
)
end


    function wheeler_diagram!(
    ax1::Axis,
    ax2::Axis,
    header::Header,
    data::DataSlice;
    smooth_size=(3,11),
    colormap=Reverse(:curl),
    range=(-100.0u"m/Myr", 100.0u"m/Myr"),
    facies_colors=nothing,
    nondep_eps=0.0u"m"
)

    linkyaxes!(ax1, ax2)

    sa = sediment_accumulation!(
        ax1, header, data;
        smooth_size=smooth_size,
        colormap=colormap,
        range=range
    )

    ft = dominant_facies!(
        ax2, header, data;
        smooth_size=smooth_size,
        facies_colors=facies_colors,
        nondep_eps=nondep_eps
    )

    ax2.ylabel = ""
    return sa, ft
end


function wheeler_diagram(header::Header, data::DataSlice; kwargs...)
    fig = Figure(size=(1000,600))
    ax1 = Axis(fig[2,1])
    ax2 = Axis(fig[2,2])
    sa, ft = wheeler_diagram!(ax1, ax2, header, data; kwargs...)
    Colorbar(fig[1,1], sa; vertical=false, label="sediment accumulation [m/Myr]")
    Colorbar(fig[1,2], ft; vertical=false, label="dominant facies")
    fig
end

end
