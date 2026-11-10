# ~/~ begin <<docs/src/visualization.md#ext/StratigraphicColumn.jl>>[init]
module StratigraphicColumn

using Makie
using Unitful

import CarboKitten.Visualization: stratigraphic_column!
using CarboKitten.Export: age_depth_model
using CarboKitten.Output.Storage: Header, DataColumn, Pack, PackColumn,
    stratigraphic_column, sediment_thickness


function scdata(p::PackColumn)
    sc = stratigraphic_column(p) / u"m"

    colormax(d) = getindex.(argmax(d; dims=1)[1, :], 1)
    adm = age_depth_model(sediment_thickness(p))

    return (ys_low=adm[1:end-1] / u"m", ys_high=adm[2:end] / u"m", facies=colormax(sc)[1:end-1])
end


function stratigraphic_column!(ax::Axis, p::PackColumn; color=Makie.wong_colors())
    (ys_low, ys_high, facies) = scdata(p)
    hspan!(ax, ys_low, ys_high; color=color[facies])
end

function stratigraphic_column!(ax::Axis, p::Observable{<:PackColumn}; color=Makie.wong_colors())
    _scdata = lift(scdata, p)
    _ys_low = lift(d -> d.ys_low, _scdata)
    _ys_high = lift(d -> d.ys_high, _scdata)
    _color = lift(d -> color[d.facies], _scdata)
    hspan!(ax, _ys_low, _ys_high; color=_color)
end

# Methods accepting a `Header`, `DataColumn` pair, kept for backward compatibility.
stratigraphic_column!(ax::Axis, header::Header, data::DataColumn; kwargs...) =
    stratigraphic_column!(ax, Pack(header, data); kwargs...)
stratigraphic_column!(ax::Axis, header::Header, data::Observable{DataColumn}; kwargs...) =
    stratigraphic_column!(ax, lift(d -> Pack(header, d), data); kwargs...)

end
# ~/~ end
