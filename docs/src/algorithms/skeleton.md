## Skeleton

The skeleton of a bitmap reduces the dimension of its structures to line segments. Our input is a `Matrix{Bool}`, where we want to plot lines over the regions where the input is `True`. With our application to visualizing hiatuses in mind, our task is made simpler. We are limited to a reduction from two dimensions down to line segments in a fixed axis.

Suppose we only needed to operate on individual columns: scanning through time we are then tasked to find consecutive regions where our input is true. This is done by the [Range Finder](range_finder.md) iterator.
To learn whether we need to plot a line, we need to find where the found ranges in neighbouring columns overlap. That is supposing we are plotting for a slice with a fixed `y` coordinate, given a column with location `x`, we use `find_ranges` to generate our vertices `(x, z1:z2)`.

``` {.julia #skeleton}
const Vertex = Tuple{Int, UnitRange{Int}}
```

We do this for all columns to find all vertices. To see if two vertices are connected by an edge, there are two criteria:

1. The vertices need to be adjacent in the `x` direction. We don't need to make this selection if we iterate over adjacent pairs directly:

   ``` {.julia #skeleton}
   pairs(it) = zip(it, drop(it, 1))
   ```

2. Their ranges should have a non-zero intersection.

   ``` {.julia #skeleton}
   edge(a::Vertex, b::Vertex) = isempty(a[2] ∩ b[2]) ? nothing : (a[1], b[1])
   ```

For two adjacent columns we need to cross-check all their ranges.

``` {.julia #skeleton}
edges_between(a, b) = filter(!isnothing, imap(splat(edge), product(a, b)))
```

To get the `z` coordinate of a vertex, we need take the average.

``` {.julia #skeleton}
middle(a::UnitRange{Int}) = (a.start + a.stop) ÷ 2
```

All this, taken together, lets us write the `skeleton` function in a quite terse, but admittedly slightly unreadable manner. We have one iterator `find_ranges` that we map over each row of the input, then a second iterator `enumerate_seq` (see [Nested Sequence Enumeration](enumerate_seq.md)) that tags all found vertices with an integer id.

``` {.julia file=src/Algorithms/Skeleton.jl}
module Skeleton

using .Iterators: filter, map as imap, product, flatten, drop
using ..Algorithms: enumerate_seq, find_ranges

<<skeleton>>

"""
    skeleton(bitmap::AbstractMatrix{Bool})

Computes the skeleton of a bitmap, i.e. reduces features with some thickness to
a set of line segments. This function is designed with stratigraphic application
in mind: we scan each row in the bitmap for connected regions, then link neighbouring
regions when they overlap. The result is a graph that represents hiatus in the sediment
accumulation.

Returns a tuple of `vertices` and `edges`, where `vertices` is a vector of 2-tuples and
`edges` is a nx2 matrix of indices into the `vertices`.
"""
function skeleton(bitmap::AbstractMatrix{Bool}; minwidth=10)
    vertex_rows = (filter(r->length(r)>=minwidth, find_ranges(row)) for row in eachrow(bitmap))
    edges::Vector{Tuple{Int,Int}} =
        map(splat(edges_between), pairs(enumerate_seq(vertex_rows))) |> flatten |> collect
    vertices::Vector{Tuple{Int,Int}} =
        (((i, middle(v)) for v in vs) for (i, vs) in enumerate(vertex_rows)) |> flatten |> collect
    return vertices, reshape(reinterpret(Int, edges), (2,:))'
end

export skeleton

end
```

## Tests

The skeleton is a graph, so we can test its topology. For a graph with $V$ vertices, $E$ edges and $C$ connected components, the number of independent cycles (loops) is $E - V + C$. In other words, $V - E$ equals the number of components minus the number of loops.

- A disc reduces to a tree: one component, no loops, so $V - E = 1$.
- An annulus reduces to a single loop: one component, one loop, so $V - E = 0$.
- Two separate discs give two trees: $V - E = 2$.

Counting $V - E$ alone doesn't tell these cases apart in general, so we also count the connected components with a small union-find.

``` {.julia file=test/Algorithms/SkeletonSpec.jl}
module SkeletonSpec
using Test
using CarboKitten.Algorithms: skeleton

centered_grid(n) = ((i - (n + 1) / 2, j - (n + 1) / 2) for i in 1:n, j in 1:n)
disc(n, r) = [hypot(x, y) <= r for (x, y) in centered_grid(n)]
annulus(n, r1, r2) = [r1 <= hypot(x, y) <= r2 for (x, y) in centered_grid(n)]

function n_components(n_vertices, edges)
    parent = collect(1:n_vertices)
    root(i) = parent[i] == i ? i : (parent[i] = root(parent[i]))
    for (a, b) in eachrow(edges)
        parent[root(a)] = root(b)
    end
    return count(i -> root(i) == i, 1:n_vertices)
end

"""
Returns the number of connected components and the number of independent
loops in the skeleton of `bitmap`.
"""
function topology(bitmap; minwidth=1)
    vertices, edges = skeleton(bitmap; minwidth=minwidth)
    v = length(vertices)
    e = size(edges, 1)
    c = n_components(v, edges)
    return (components=c, loops=e - v + c)
end

@testset "Algorithms/Skeleton" begin
    @testset "graph structure" begin
        vertices, edges = skeleton(annulus(41, 8, 18); minwidth=1)
        @test all(1 .<= edges .<= length(vertices))
        # edges only connect vertices in adjacent rows
        @test all(vertices[b][1] - vertices[a][1] == 1 for (a, b) in eachrow(edges))
    end

    @testset "topology" begin
        @test topology(disc(41, 15)) == (components=1, loops=0)
        @test topology(annulus(41, 8, 18)) == (components=1, loops=1)
        @test topology(annulus(41, 15, 17)) == (components=1, loops=1)
        @test topology(hcat(disc(21, 8), disc(21, 8))) == (components=2, loops=0)
    end

    @testset "minwidth" begin
        # The walls of the annulus are about 11 cells wide. A larger `minwidth`
        # removes the sides, leaving the top and bottom caps as two trees.
        @test topology(annulus(41, 8, 18); minwidth=12) == (components=2, loops=0)
        @test isempty(skeleton(annulus(41, 8, 18); minwidth=40)[1])
    end

    @testset "empty input" begin
        vertices, edges = skeleton(falses(5, 5); minwidth=1)
        @test isempty(vertices)
        @test size(edges) == (0, 2)
    end
end
end
```
