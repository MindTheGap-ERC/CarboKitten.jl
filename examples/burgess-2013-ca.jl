# ~/~ begin <<docs/src/carbocat-ca.md#examples/burgess-2013-ca.jl>>[init]
using CarboKitten
using CarboKitten.Burgess2013.CA
using CarboKitten.Stencil: Reflected
using CarboKitten.Utility
using Plots

function main()
    plotlyjs()
    l = @layout([a b c d; e f g h])

    init = rand(0:3, 50, 50)
    result = Iterators.take(CA.run(Reflected{2}, init, 3), 8)

    plot((heatmap(r, colorbar=:none, 
                     aspect_ratio=1, 
                     xlims=(0, 50), 
                     ylims=(0, 50)) for r in result)...,
        layout=(2, 4),
        size=(800, 400))
    savefig("docs/src/fig/b13-first-8-iterations.html")
end


# plot("burgess-fig3.svg")

function plot_long_times(output::String)
    init = rand(0:3, 50, 50)
    result = select(CA.run(Reflected{2}, init, 3), [10, 100, 1000])
end

# plot_long_times("burgess-long.svg")

main()
# ~/~ end