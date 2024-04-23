# ~/~ begin <<docs/src/stencils.md#examples/ca/life.jl>>[init]
#| creates: docs/src/fig/life.gif
#| requires: src/Stencil.jl
#| collect: figures

module Life
    using CarboKitten.Stencil
    using GLMakie
    using .Iterators: take

    "x is a 3x3 region around the cell at x[2,2]."
    rules(x) = let c = x[2, 2], s = sum(x) - c
        c && s == 2 || s == 3
    end

    function game_of_life(w, h)
        y1 = rand(Bool, (w, h))
        y2 = Array{Bool}(undef, w, h)

        op = stencil(Bool, Periodic{2}, (3, 3), rules)
        Channel() do ch
            put!(ch, y1)
            while true
                op(y1, y2)
                (y1, y2) = (y2, y1)
                put!(ch, y1)
            end
        end
    end

    function plot()
        life = take(game_of_life(50, 50), 150)
        fig = Figure()
        ax = Axis(fig[1,1], aspect=1)
        record(fig, "docs/src/fig/life.gif", life; framerate=10) do frame
            heatmap!(ax, frame; colormap=:Blues)
        end
    end
end

Life.plot()
# ~/~ end