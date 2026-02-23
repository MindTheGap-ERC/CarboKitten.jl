# ~/~ begin <<docs/src/components/production.md#examples/production/pelagic.jl>>[init]
module PelagicProductionPlot

using CarboKitten
using CarboKitten.Components.Production: benthic_production, pelagic_production
using CarboKitten.ALCAP.Example: INPUT
using CairoMakie

function main()
    water_depth = (0.01:1.0:70.0)u"m"
    fig = Figure()
    facies = [f for f in INPUT.facies if f.active]

    ax_benthic = Axis(fig[1, 1], yreversed=true, title="benthic production", ylabel="depth [m]")
    for f in facies
        p = water_depth .|> (d -> benthic_production(INPUT.insolation, f, d))
        lines!(ax_benthic, p |> in_units_of(u"m/Myr"),
            water_depth |> in_units_of(u"m"), label = string(f.name))
    end

    ax_pelagic = Axis(fig[1, 2], yreversed=true, title="pelagic production")
    for f in facies
        p = water_depth .|> (d -> pelagic_production(INPUT.insolation, f, d))
        lines!(ax_pelagic, p |> in_units_of(u"m^2/Myr"),
            water_depth |> in_units_of(u"m"), label = string(f.name))
    end
    fig[1, 3] = Legend(fig, ax_benthic, "Facies")
    fig
end

end
# ~/~ end
