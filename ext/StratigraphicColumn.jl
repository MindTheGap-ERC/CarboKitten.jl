# ~/~ begin <<docs/src/visualization.md#ext/StratigraphicColumn.jl>>[init]
module StratigraphicColumn

using Makie
using Unitful

import CarboKitten                      
import CarboKitten.Visualization: stratigraphic_column!, stratigraphic_column_layers!, facies_colormap

using CarboKitten.Export: datacolumn_from_state
using CarboKitten.Export: Header, DataColumn, stratigraphic_column, age_depth_model

export stratigraphic_column_layers!

function CarboKitten.Visualization.stratigraphic_column_layers!(
    ax::Makie.Axis,
    state;
    x_index::Int = 1,
    y_index::Int = 1,
    facies_colors = nothing
)
    comp = state.compaction

    n = comp.n_layers[x_index, y_index]
    n == 0 && return

n_facies = length(facies_colors)         

cmap = facies_colormap(
    n_facies;
    facies_colors = facies_colors,
    include_nodeposit = false
)

colors = Makie.to_color.(cmap.colors)

    z = 0.0

    for k in reverse(1:n)
        f = comp.layer_facies[x_index, y_index, k]
        h = comp.layer_thickness[x_index, y_index, k]

        if f != 0   # skip NoDeposit
    hspan!(ax, [z], [z + h];
        color = colors[f]
    )
end

        z += h
    end

    ax.ylabel = "Depth (m)"
    ax.xlabel = "Facies"
    ax.yreversed = true
end

function stratigraphic_column!(ax::Makie.Axis, header::Header, state;
    x_index = 1,
    y_index = 1,
    write_interval = 1,
    facies_colors = nothing
)
    col = datacolumn_from_state(
        state;
        x_index = x_index,
        y_index = y_index,
        write_interval = write_interval
    )

    stratigraphic_column!(ax, header, col; facies_colors = facies_colors)
end
end
