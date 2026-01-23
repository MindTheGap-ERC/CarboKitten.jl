# ~/~ begin <<docs/src/algorithms/skeleton.md#src/Algorithms/Skeleton.jl>>[init]
module Skeleton

using .Iterators: filter, map as imap, product, flatten, drop
using ..Algorithms: enumerate_seq, find_ranges

# ~/~ begin <<docs/src/algorithms/skeleton.md#skeleton>>[init]
const Vertex = Tuple{Int, UnitRange{Int}}
# ~/~ end
# ~/~ begin <<docs/src/algorithms/skeleton.md#skeleton>>[1]
pairs(it) = zip(it, drop(it, 1))
# ~/~ end
# ~/~ begin <<docs/src/algorithms/skeleton.md#skeleton>>[2]
edge(a::Vertex, b::Vertex) = isempty(a[2] ∩ b[2]) ? nothing : (a[1], b[1])
# ~/~ end
# ~/~ begin <<docs/src/algorithms/skeleton.md#skeleton>>[3]
edges_between(a, b) = filter(!isnothing, imap(splat(edge), product(a, b)))
# ~/~ end
# ~/~ begin <<docs/src/algorithms/skeleton.md#skeleton>>[4]
middle(a::UnitRange{Int}) = (a.start + a.stop) ÷ 2
# ~/~ end

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
# ~/~ end
