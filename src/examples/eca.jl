# ~/~ begin <<docs/src/stencils.md#src/examples/eca.jl>>[init]
using CarboKitten.Stencil
using GnuplotLite

rule(i::Int) = function (foo::AbstractVector{T}) where T <: Integer
    d = foo[1]*4 + foo[2]*2 + foo[3]
    i & (1 << d) == 0 ? 0 : 1
end

function eca(r::Int, n::Int, iter::Int)
    y = Array{Int}(undef, n, iter)
    y[:, 1] = rand(0:1, n)
    stencil_op = stencil(Int, Periodic{1}, (3,), rule(r))
    for i in 2:iter
        stencil_op(view(y, :, i-1), view(y, :, i))
    end
    y
end

# ~/~ begin <<docs/src/stencils.md#eca-plot>>[init]
function plot_eca(output::String, r::Int, n::Int, iter::Int)
    plot_width = 700
    plot_height = plot_width * iter ÷ n + 100
    gnuplot() do g
        g |>
            send("set term svg size $(plot_width), $(plot_height)") |>
            send("set output '$(output)'") |>
            send("data" => (x=0:n-1, y=0:iter-1, z=eca(r, n, iter)')) |>
            send("set title 'rule $(r)'") |>
            send("set xrange [0:$(n)]; set yrange [$(iter):0] reverse") |>
            send("set xlabel 'space'") |>
            send("set ylabel 'iterations'") |>
            send("set size ratio -1") |>
            send("unset colorbox; set palette gray") |>
            send("plot \$data nonuniform matrix u (\$1+0.5):(\$2+0.5):3 with image")
    end
end
# ~/~ end
# ~/~ end