# ~/~ begin <<docs/src/algorithms/multi-threading.md#test/Interfaces/ChunksSpec.jl>>[init]
@testset "CarboKitten.Interfaces.Chunks" begin
    using CarboKitten.Interfaces.Chunks

    @test Chunk(:, 4) isa Chunk{2}
    @test !(Chunk(1:25, 50:75) isa TotalChunk)
    @test Chunk(:, :) isa TotalChunk
end
# ~/~ end
