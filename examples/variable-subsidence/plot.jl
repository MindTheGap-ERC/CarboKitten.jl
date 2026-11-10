# ~/~ begin <<docs/src/cases/variable-subsidence.md#examples/variable-subsidence/plot.jl>>[init]

module PlotVariableSubsidence

using CairoMakie
using CarboKitten
using CarboKitten.Output: total_subsidence
using CarboKitten.Visualization: sediment_profile!

const FILENAME = "data/output/variable-subsidence.h5"

function main()
    bundle = load(FILENAME)
    fig = Figure(size=(1200, 800))

    profiles = [:edge => "edge", :flank => "flank", :centre => "depocentre"]
    slices = [load_slice(bundle, k) for (k, _) in profiles]
    hdr = slices[1].header

    ax0 = Axis(fig[1, 1]; title="total subsidence", xlabel="x [km]", ylabel="y [km]",
               aspect=DataAspect())
    hm = heatmap!(ax0, hdr.axes.x |> in_units_of(u"km"), hdr.axes.y |> in_units_of(u"km"),
                  total_subsidence(hdr) |> in_units_of(u"m"); colormap=:lajolla)
    for (p, _) in zip(slices, profiles)
        y = p.header.axes.y[p.data.slice[2]] |> in_units_of(u"km")
        hlines!(ax0, y; color=:black, linestyle=:dash)
    end
    Colorbar(fig[1, 2], hm; label="subsidence [m]")

    positions = [fig[1, 3], fig[2, 1:2], fig[2, 3]]
    for (pos, p, (_, title)) in zip(positions, slices, profiles)
        ax = Axis(pos)
        sediment_profile!(ax, p)
        ax.title = title
    end

    close(bundle)
    save("docs/src/_fig/variable-subsidence.png", fig)
end

end

PlotVariableSubsidence.main()
# ~/~ end
