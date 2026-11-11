# ~/~ begin <<docs/src/production/pelagic.md#examples/production/pelagic.jl>>[init]
module PelagicProductionPlot

using CarboKitten
using CarboKitten.Production
using CarboKitten.Components: Production as ProductionComponent
using CarboKitten.Components.Common: record_type
using CairoMakie

const Input = record_type(:Input, ProductionComponent)

function main()
    water_depth = (0.01:0.1:50.0)u"m"
    fig = Figure(size=(600,600))
    input = Input(
        box = CarboKitten.Box{Coast}(grid_size=(1, 1), phys_scale=1.0u"m"),
        time = TimeProperties(Δt=1.0u"kyr", steps=1))

    ax = Axis(fig[1, 1], yreversed=true, ylabel="depth [m]", xlabel="production [m/Myr]")
    for (k, prod) in pairs(Production.EXAMPLE)
        f = production_profile(input, prod)
        p = water_depth .|> (w -> f(input.time.t0, w))
        lines!(ax, p |> in_units_of(u"m/Myr"),
            water_depth |> in_units_of(u"m"), label = string(k))
    end
    fig[1, 2] = Legend(fig, ax, "Facies")
    fig
end

end
# ~/~ end
