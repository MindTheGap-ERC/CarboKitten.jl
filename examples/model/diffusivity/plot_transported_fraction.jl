# ~/~ begin <<docs/src/active-layer-transport.md#examples/model/diffusivity/plot_transported_fraction.jl>>[init]
module Transported_fraction
using GLMakie
using CarboKitten.Export: read_slice
using CarboKitten.Visualization: profile_plot!

function main()
    (header, slice) = read_slice("data/output/diffusivity-example.h5", :profile)
    fig = Figure()
    ax = Axis(fig[1, 1])

    x = header.axes.x
#    t = header.axes.t

    plot = profile_plot!(x -> x[4]/sum(x), ax, header, slice; colorrange=(0, 0.1))
    Colorbar(fig[1, 2], plot; label=L"f_1transported / f_{total}")

    save("docs/src/_fig/transported_fraction4.png", fig)
    fig
end
end

Transported_fraction.main()
# ~/~ end
