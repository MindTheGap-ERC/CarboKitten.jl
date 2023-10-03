# ~/~ begin <<docs/src/carbocat-ca.md#examples/ca/burgess-2013.jl>>[init]
module Script
    using .Iterators: flatten
    using CarboKitten
    using CarboKitten.Burgess2013
    using CarboKitten.Stencil: Reflected
    using CarboKitten.Utility
    using GLMakie

    function main()
        init = rand(0:3, 50, 50)
        ca = run_ca(Reflected{2}, MODEL1, init, 3)

        fig = Figure(resolution=(1000, 500))
        axis_indices = flatten(eachrow(CartesianIndices((2, 4))))
        for (i, st) in zip(axis_indices, ca)
            ax = Axis(fig[Tuple(i)...], aspect=AxisAspect(1))
            heatmap!(ax, st)
        end
        save("docs/src/fig/b13-fig3.png", fig)
    end
end

Script.main()
# ~/~ end