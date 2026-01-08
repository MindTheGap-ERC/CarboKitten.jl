# ~/~ begin <<docs/src/algorithms/enumerate_seq.md#src/Algorithms/EnumerateSeq.jl>>[init]
module EnumerateSeq

struct TagVectors{T}
    vectors::T
end

function Base.iterate(tv::TagVectors{T}) where T
    tag = 1
    x = iterate(tv.vectors)
    isnothing(x) && return nothing
    (v_, nit) = x
    v = collect(v_)
    n = length(v)
    return zip(tag:tag+n-1, v) |> collect, (tag+n, nit)
end

function Base.iterate(tv::TagVectors{T}, st::Tuple{Int,U}) where {T, U}
    (tag, it) = st
    isnothing(it) && return nothing
    x = iterate(tv.vectors, it)
    isnothing(x) && return nothing
    (v_, nit) = x
    v = collect(v_)
    n = length(v)
    return zip(tag:tag+n-1, v) |> collect, (tag+n, nit)
end

Base.IteratorSize(::Type{TagVectors{T}}) where T = Base.IteratorSize(T)
Base.size(r::TagVectors{T}) where T = size(r.vectors)
Base.length(r::TagVectors{T}) where T = length(r.vectors)
Base.eltype(::Type{TagVectors{T}}) where T = Any
# This should be:
# Iterators.Zip{Tuple{UnitRange{Int},Vector{eltype(eltype(T))}}}
# But that doesn't work

"""
    enumerate_seq(s)

Enumerates an iterator of sequences, such that the following equivalence
holds,

    enumerate(flatten(s)) == flatten(enumerate_seq(s))

while keeping the underlying structure intact.
"""
enumerate_seq(s::T) where T = TagVectors{T}(s)

export enumerate_seq

end
# ~/~ end
