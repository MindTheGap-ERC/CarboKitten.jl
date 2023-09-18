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

## Howto run
We start with randomized initial conditions on a 50x50 grid.

``` {.julia}
init = rand(0:3, 50, 50)
```

Then we run the cellular automaton for, in this case eight generations. The `CA.run` function returns an iterator of 50x50 maps. That means that in principle we can extract an infinity of iterations, but in this cane we only `take` eight.

``` {.julia}
result = Iterators.take(CA.run(Reflected{2}, init, 3), 8)
```

In Julia we may plot those as follows

``` {.julia}
using Plots
# plotly()  # sets back-end; plotly gives me the best results
plot((heatmap(r, colorbar=:none) for r in result)..., layout=(2, 4))
```

What this says is: create a `heatmap` for each of our eight results, then expand those into a function call to `plot` (as in `plot(hm1, hm2, ..., hm8, layout=(2, 4))`).

<details><summary>Plotting code</summary>

``` {.julia file=examples/burgess-2013-ca.jl}
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
```

</details>

