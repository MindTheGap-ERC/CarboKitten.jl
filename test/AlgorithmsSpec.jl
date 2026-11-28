# ~/~ begin <<docs/src/algorithms/overview.md#test/AlgorithmsSpec.jl>>[init]
module AlgorithmsSpec
using Test
using CarboKitten.Algorithms: enumerate_seq

@testset "CarboKitten.Algorithms" begin
    # ~/~ begin <<docs/src/algorithms/enumerate_seq.md#algorithms-spec>>[init]
    a = [[:a, :b], [:c]]
    @test collect(enumerate_seq(a)) == [[(1, :a), (2, :b)], [(3, :c)]]
    # ~/~ end
end
end
# ~/~ end
