# ~/~ begin <<docs/src/visualization.md#ext/WheelerDiagram.jl>>[init]
module WheelerDiagram

import CarboKitten.Visualization: wheeler_diagram!, wheeler_diagram
using CarboKitten.Export: Header, Data, DataSlice, read_data, read_slice
using CarboKitten.Utility: in_units_of
using Makie
using Unitful

elevation(h::Header, d::DataSlice) =
    let bl = h.bedrock_elevation[d.slice..., na],
        sr = h.axes.t[end] * h.subsidence_rate

        cat(bl, bl .+ d.sediment_elevation; dims=2) .- sr
    end

function wheeler_diagram!(ax::Axis, header::Header, data_slice::DataSlice)
    colormax(d) = getindex.(argmax(d; dims=1)[1, :, :], 1)
    magnitude = sum(data_slice.deposition .- data_slice.disintegration; dims=1)[1, :, :]
    dominant_facies = colormax(data_slice.deposition)

    ξ = elevation(header, data_slice)  # |> in_units_of(u"m")
    water_depth = ξ .- (header.subsidence_rate.*(header.axes.t.-header.axes.t[end]).+header.sea_level)[na, :]
    exposed = water_depth .< 0.0u"m"

    heatmap!(ax, dominant_facies)
    contourf!(ax, exposed)
end

end
# ~/~ end
