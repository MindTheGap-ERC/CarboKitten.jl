# ~/~ begin <<docs/src/carbocat-ca.md#examples/ca/long-term.jl>>[init]
using CarboKitten
using CarboKitten.Burgess2013
using CarboKitten.Stencil
using CarboKitten.Utility
using GLMakie

function main()
    init = rand(0:3, 50, 50)
    result = select(run_ca(Periodic{2}, MODEL1, init, 3), [10, 100, 10000])

    fig = Figure(resolution=(1000, 333))
    for (i, st) in enumerate(result)
        ax = Axis(fig[1, i], aspect=AxisAspect(1))
        heatmap!(ax, st)
    end
    save("docs/src/fig/b13-long-term.png", fig)
end

main()
# ~/~ end