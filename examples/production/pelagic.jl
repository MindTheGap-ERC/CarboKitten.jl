# ~/~ begin <<docs/src/components/production.md#examples/production/pelagic.jl>>[init]
module PelagicProductionPlot

using CarboKitten
using CarboKitten.Production
using CairoMakie

@kwdef struct Input <: CarboKitten.AbstractInput
    insolation = 400.0u"W/m^2"
end

function main()
    water_depth = (0.01:0.1:50.0)u"m"
    fig = Figure(size=(600,600))
    input = Input()

    ax = Axis(fig[1, 1], yreversed=true, ylabel="depth [m]", xlabel="production [m/Myr]")
    for (k, prod) in pairs(Production.EXAMPLE)
        f = production_profile(input, prod)
        p = water_depth .|> (w -> f(input.insolation, w))
        lines!(ax, p |> in_units_of(u"m/Myr"),
            water_depth |> in_units_of(u"m"), label = string(k))
    end
    fig[1, 2] = Legend(fig, ax, "Facies")
    fig
end

end
# ~/~ end
