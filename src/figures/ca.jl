# ~/~ begin <<docs/src/carbocat-ca.md#src/figures/ca.jl>>[init]
using MindTheGap.Burgess2013.CA
using MindTheGap.Stencil: Reflected
using MindTheGap.Utility
using GnuplotLite

function plot_array(w::Int, h::Int, msgs; xtics="set xtics", ytics="set ytics")
    ch = Channel{String}() do ch
        g = Gnuplot(ch)

        top_margin = 0.002
        bottom_margin = 0.1
        left_margin = 0.05
        right_margin = 0.05
        inner_margin = 0.02
        vert_inner_margin = 0.025

        # h * ph + (h-1) * inner_margin + top_margin + bottom_margin = 1
        plot_height = (1.0 - (h-1)*inner_margin - top_margin - bottom_margin) / h
        plot_width = (1.0 - (w-1)*vert_inner_margin - left_margin - right_margin) / w

        g |> send("set multiplot; unset xtics")
        for (i, msg) in enumerate(msgs)
            grid_x = (i - 1) % w
            grid_y = (i - 1) ÷ w
            if grid_x == 0
                tmargin = 1.0 - top_margin - (plot_height + vert_inner_margin) * grid_y
                bmargin = tmargin - plot_height
                g |> send("set tmargin at screen $(tmargin); set bmargin at screen $(bmargin)")
                g |> send(ytics)
            end
            if grid_x == 1
                g |> send("unset ytics")
            end
            if grid_y == h-1 && grid_x == 0
                g |> send(xtics)
            end
            lmargin = left_margin + (plot_width + inner_margin) * grid_x
            rmargin = lmargin + plot_width
            g |> send("set lmargin at screen $(lmargin); set rmargin at screen $(rmargin)")
            g |> msg
        end
        g |> send("unset multiplot")
    end

    send(join(ch, "\n"))
end

function plot(output::String)
    init = rand(0:3, 50, 50)
    result = Iterators.take(CA.run(Reflected{2}, init, 3), 8)

    gnuplot() do g
        g |>
            send("set term svg size 900, 440") |>
            send("set output '$(output)'") |>
            send("load 'data/blue-to-red.pal'") |>
            send("set size square") |>
            send("set xrange [0:50]; set yrange [0:50]") |>
            send("unset colorbox") |>
            plot_array(4, 2, (send("data" => r) *
                              send("plot \$data matrix u (\$1+0.5):(\$2+0.5):3 t'' w image pixels")
                              for r in result);
                       xtics="set xtics 10, 10, 50",
                       ytics="set ytics 10, 10, 50")
    end
end

function plot_long_times(output::String)
    init = rand(0:3, 50, 50)
    result = select(CA.run(Reflected{2}, init, 3), [10, 100, 1000])

    gnuplot() do g
        g |>
            send("set term svg size 900, 240") |>
            send("set output '$(output)'") |>
            send("load 'data/blue-to-red.pal'") |>
            send("set size square") |>
            send("set xrange [0:50]; set yrange [0:50]") |>
            send("unset colorbox") |>
            plot_array(3, 1, (send("data" => r) *
                              send("plot \$data matrix u (\$1+0.5):(\$2+0.5):3 t'' w image pixels")
                              for r in result);
                       xtics="set xtics 10, 10, 50",
                       ytics="set ytics 10, 10, 50")
    end
end
# ~/~ end