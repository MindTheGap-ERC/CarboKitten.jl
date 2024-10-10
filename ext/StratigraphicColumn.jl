# ~/~ begin <<docs/src/visualization.md#ext/StratigraphicColumn.jl>>[init]
module StratigraphicColumn

using Makie
using Unitful

import CarboKitten.Visualization: stratigraphic_column!
using CarboKitten.Export: Header, DataColumn, stratigraphic_column, age_depth_model


function scdata(header::Header, data::DataColumn)
    n_facies = size(data.production)[1]
    n_times = length(header.axes.t) - 1
    sc = zeros(Float64, n_facies, n_times)
    for f = 1:n_facies
        sc[f, :] = stratigraphic_column(header, data, f) / u"m"
    end

    colormax(d) = getindex.(argmax(d; dims=1)[1, :], 1)
    adm = age_depth_model(data.sediment_elevation)

    return (ys_low=adm[1:end-1] / u"m", ys_high=adm[2:end] / u"m", facies=colormax(sc)[1:end-1])
end


function stratigraphic_column!(ax::Axis, header::Header, data::DataColumn; color=Makie.wong_colors())
    (ys_low, ys_high, facies) = scdata(header, data)
    hspan!(ax, ys_low, ys_high; color=color[facies])
end

function stratigraphic_column!(ax::Axis, header::Header, data::Observable{DataColumn}; color=Makie.wong_colors())
    _scdata = lift(d -> scdata(header, d), data)
    _ys_low = lift(d -> d.ys_low, _scdata)
    _ys_high = lift(d -> d.ys_high, _scdata)
    _color = lift(d -> color[d.facies], _scdata)
    hspan!(ax, _ys_low, _ys_high; color=_color)
end

end
# ~/~ end