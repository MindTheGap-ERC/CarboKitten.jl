---
title: Species Habitation
subtitle: a cellular automaton
---

## Cellular Automaton
The paper talks about cycling the order of preference for occupying an empty cell at each iteration. This means that the rules change slighly every iteration.

``` {.julia #cycle-permutation}
cycle_permutation(n_species::Int) =
    (circshift(1:n_species, x) for x in Iterators.countfrom(0))
```

The `stencil` function has an `args...` variadic arguments that are forwarded to the given rule. This means we can create a `rules` function that we pass the preference order as a second argument.

``` {.julia #burgess2013-rules}
function rules(neighbourhood::Matrix{Int}, order::Vector{Int})
    cell_species = neighbourhood[3, 3]
    neighbour_count(species) = sum(neighbourhood .== species)
    if cell_species == 0
        for species in order
            n = neighbour_count(species)
            if 6 <= n && n <= 10
                return species
            end
        end
        0
    else
        n = neighbour_count(cell_species) - 1
        (4 <= n && n <= 10 ? cell_species : 0)
    end
end
```

This function is not yet adaptible to the given rule set. Such a modification is not so hard to make. 

The paper talks about a 50x50 grid initialized with uniform random values.

``` {.julia file=src/Burgess2013/CA.jl}
module CA

using ...Stencil

<<cycle-permutation>>
<<burgess2013-rules>>

function run(::Type{B}, init::Matrix{Int}, n_species::Int) where {B <: Boundary{2}}
    Channel{Matrix{Int}}() do ch
        target = Matrix{Int}(undef, size(init))
        put!(ch, init)
        stencil_op = stencil(Int, B, (5, 5), rules)
        for perm in cycle_permutation(n_species)
            stencil_op(init, target, perm)
            init, target = target, init
            put!(ch, init)
        end
    end
end

end
```

First, let us reproduce Figure 3 in Burgess 2013.

![First 8 generations](burgess-fig3.svg)

By eye comparison seems to indicate that this CA is working the same. I'm curious to the behaviour after more iterations. Let's try 10, 100, 1000 and so on.

![Assymptotic behaviour](burgess-long.svg)

The little qualitative change between 100 and 1000 iterations would indicate that this CA remains "interesting" for a long time.

On my laptop I can run about 150 iterations per second with current code. When using periodic boundaries, I get to 1500 iterations per second, which is peculiar. A lot can still be optimized.

<details><summary>Plotting code</summary>

```@example
using CarboKitten
using CarboKitten.Burgess2013.CA
using CarboKitten.Stencil: Reflected
using CarboKitten.Utility
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

plot("burgess-fig3.svg")

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

plot_long_times("burgess-long.svg")
```

</details>

