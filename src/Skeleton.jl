# ~/~ begin <<docs/src/visualization.md#src/Skeleton.jl>>[init]
module Skeleton

using .Iterators: filter, map as imap, product, flatten, drop
using ..Utility: enumerate_seq, find_ranges

const Vertex = Tuple{Int, UnitRange{Int}}

pairs(it) = zip(it, drop(it, 1))
edge(a::Vertex, b::Vertex) = isempty(a[2] ∩ b[2]) ? nothing : (a[1], b[1])
edges_between(a, b) = filter(!isnothing, imap(splat(edge), product(a, b)))
middle(a::UnitRange{Int}) = (a.start + a.stop) ÷ 2


function skeleton(bitmap::AbstractMatrix{Bool}; minwidth=10)
    vertex_rows = (filter(r->length(r)>=minwidth, find_ranges(row)) for row in eachrow(bitmap))
    edges::Vector{Tuple{Int,Int}} = collect(flatten(map(splat(edges_between), pairs(enumerate_seq(vertex_rows)))))
    vertices::Vector{Tuple{Int,Int}} = collect(flatten(((i, middle(v)) for v in vs) for (i, vs) in enumerate(vertex_rows)))
    return vertices, reshape(reinterpret(Int, edges), (2,:))'
end

end
# ~/~ end