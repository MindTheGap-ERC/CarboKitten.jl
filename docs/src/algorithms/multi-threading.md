Multi-threading
===============

Currently CarboKitten is single-threaded, but it can be made to run multi-threaded. Several strategies:

Higher Level Abstractions
-------------------------

Use `Threads.@threads :static` for-loops with automated chunking, or `@floops` or some other abstraction everywhere. This adds overhead of starting multiple threads in many different places. There are quite a few stencil operations in the code which then require synchronization before continuing with the next step. For instance, when we integrate the transport equation, we need to synchronize after each transport loop. If we keep on fireing new tasks every loop, the overhead of creating and synchronizing these tasks may become too large.

Manual chunking
---------------

The user may specify a list of grid chunks over which the model stepper can parallelize. We can provide a convenience routine that generates nice chunks for the available number of threads. We may use [SyncBarriers.jl](https://github.com/JuliaConcurrent/SyncBarriers.jl) to introduce barriers in the code. Each thread is responsible for a fixed area in the simulation volume. Each routine in CarboKitten should be made to receive an extra argument indicating the chunk it is supposed to work on. For example,

```julia
push_sediment(buffer::CircularBuffer, sediment::AbstractArray{F, 3}, chunk::Chunk{2})
```

where `Chunk` is something of the form:

``` {.julia file=src/Interfaces/Chunks.jl}
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
```

``` {.julia file=test/Interfaces/ChunksSpec.jl}
@testset "CarboKitten.Interfaces.Chunks" begin
    using CarboKitten.Interfaces.Chunks

    @test Chunk(:, 4) isa Chunk{2}
    @test !(Chunk(1:25, 50:75) isa TotalChunk)
    @test Chunk(:, :) isa TotalChunk
end
```

This means we can dispatch on whether we are running multi-threaded or not.

Decision
--------

Both of these approaches require some architecting. The manual approach has the advantage that we can overload essential routines to work in both chunked and unchunked versions, allowing for gradual introduction. If we start to sprinkle parallel loops around the code for the first approach, it may be harder to switch off.
