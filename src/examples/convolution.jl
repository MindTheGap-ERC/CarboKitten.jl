# ~/~ begin <<docs/src/stencils.md#src/examples/convolution.jl>>[init]
using CarboKitten.Stencil
using GnuplotLite

function plot_boundary_types(output::String)
    n = 16
    y0 = zeros(Float64, n, n)
    y0[1, 1] = 1
    y0[n, n] = 2
    x = collect(-2:0.25:2)
    k = exp.(-(x.^2 .+ x'.^2))
    k ./= sum(k)

    y_periodic = Array{Float64}(undef, n, n)
    convolution(Periodic{2}, k)(y0, y_periodic)
    y_reflected = Array{Float64}(undef, n, n)
    convolution(Reflected{2}, k)(y0, y_reflected)
    y_constant = Array{Float64}(undef, n, n)
    convolution(Constant{2, 0.1}, k)(y0, y_constant)

    gnuplot() do g
        g |>
            send("load 'data/moreland.pal'") |>
            send("periodic" => y_periodic) |>
            send("reflected" => y_reflected) |>
            send("constant" => y_constant) |>
            send("set term svg size 700, 300") |>
            send("set output '$(output)'") |>
            send("set multiplot layout 1, 3") |>
            send("set size square") |>
            send("set xrange [0:$(n)]; set yrange [0:$(n)]") |>
            send("unset colorbox") |>
            # send("set colorbox horiz user origin graph 0,screen .04 size graph 1,screen .04") |>
            send("unset xtics; unset ytics") |>
            send("set title 'periodic'") |>
            send("plot \$periodic matrix u (\$1+0.5):(\$2+0.5):3 t'' w image") |>
            send("set title 'reflected'") |>
            send("plot \$reflected matrix u (\$1+0.5):(\$2+0.5):3 t'' w image") |>
            send("set title 'constant (0.1)'") |>
            send("plot \$constant matrix u (\$1+0.5):(\$2+0.5):3 t'' w image") |>
            send("unset multiplot")
    end
end
# ~/~ end