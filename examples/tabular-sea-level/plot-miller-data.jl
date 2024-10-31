# ~/~ begin <<docs/src/cases/tabular-sea-level.md#examples/tabular-sea-level/plot-miller-data.jl>>[init]
#| creates: docs/src/_fig/miller-sea-leavel.svg
module PlotMillerData

using CarboKitten.Components.Common
using CarboKitten.DataSets: miller_2020
using CairoMakie
using DataFrames
using Unitful

function main()
    df = miller_2020()
    fig = Figure(size=(1000,300))
    ax = Axis(fig[1, 1]; xlabel="time (Ma BP)", ylabel="sealevel (m)")

    for ref in levels(df.reference)
        subset = df[df.reference .== ref,:]
        lines!(ax, subset.time |> in_units_of(u"Myr"), subset.sealevel |> in_units_of(u"m"), label=ref)
    end
    fig[1, 2] = Legend(fig, ax)

    save("docs/src/_fig/miller-sea-level.svg", fig)
    fig
end

end
# ~/~ end
