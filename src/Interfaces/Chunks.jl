# ~/~ begin <<docs/src/algorithms/multi-threading.md#src/Interfaces/Chunks.jl>>[init]
module Chunks

export Serial, Parallel, Chunk, NormalChunk, TotalChunk, normalize

const Selector = Union{Colon, UnitRange{Int64}, Int64}
const Chunk{N} = NTuple{N, Selector}
const NormalChunk{N} = NTuple{N, UnitRange{Int64}}
const TotalChunk{N} = NTuple{N, Colon}

Chunk(args::Vararg{Selector, N}) where N = Chunk{N}((args...,))

normalize(chunk::NormalChunk, args...) = chunk

function normalize(chunk::TotalChunk, a::AbstractArray)
    (_..., x, y) = size(a)
    return Chunk(1:x, 1:y)
end

struct Serial{C}
    slice::C
end

struct Parallel{C}
    slice::C
end

end
# ~/~ end
