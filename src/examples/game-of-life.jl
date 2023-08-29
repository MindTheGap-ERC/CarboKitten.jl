# ~/~ begin <<docs/src/stencils.md#src/examples/game-of-life.jl>>[init]
using CarboKitten.Stencil
using GnuplotLite

function game_of_life(w, h)
    y1 = rand(Bool, (w, h))
    y2 = Array{Bool}(undef, w, h)

    # ~/~ begin <<docs/src/stencils.md#game-of-life-rules>>[init]
    "x is a 3x3 region around the cell at x[2,2]."
    rules(x) = let c = x[2, 2], s = sum(x) - c
        c && s == 2 || s == 3
    end
    # ~/~ end

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

# ~/~ begin <<docs/src/stencils.md#life-plot>>[init]
function plot_life(output::String, w::Int, h::Int)
    (z, _) = Iterators.peel(Iterators.drop(game_of_life(w, h), 50))
    plot_width = 700
    plot_height = plot_width * h ÷ w + 100
    gnuplot() do g
        g |>
            send("set term svg size $(plot_width), $(plot_height)") |>
            send("set output '$(output)'") |>
            send("data" => Matrix{Int}(z')) |>
            send("set title 'game of life'") |>
            send("set xrange [0:$(w)]; set yrange [0:$(h)]") |>
            send("set size ratio -1") |>
            send("unset colorbox; set palette gray negative") |>
            send("plot \$data matrix u (\$1+0.5):(\$2+0.5):3 with image pixels")
    end
end
# ~/~ end
# ~/~ end