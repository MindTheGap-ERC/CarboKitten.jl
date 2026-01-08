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
