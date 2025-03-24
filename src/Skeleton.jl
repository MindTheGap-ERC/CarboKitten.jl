# ~/~ begin <<docs/src/visualization.md#src/Skeleton.jl>>[init]
module Skeleton

using .Iterators: filter, map as imap, product, flatten, drop
using ..Utility: enumerate_seq, find_ranges

const Vertex = Tuple{Int, UnitRange{Int}}

pairs(it) = zip(it, drop(it, 1))
edge(a::Vertex, b::Vertex) = isempty(a[2] ∩ b[2]) ? nothing : (a[1], b[1])
edges_between(a, b) = filter(!isnothing, imap(splat(edge), product(a, b)))
middle(a::UnitRange{Int}) = (a.start + a.stop) ÷ 2

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
    edges::Vector{Tuple{Int,Int}} = collect(flatten(map(splat(edges_between), pairs(enumerate_seq(vertex_rows)))))
    vertices::Vector{Tuple{Int,Int}} = collect(flatten(((i, middle(v)) for v in vs) for (i, vs) in enumerate(vertex_rows)))
    return vertices, reshape(reinterpret(Int, edges), (2,:))'
end

end
# ~/~ end
