# ~/~ begin <<docs/src/algorithms/overview.md#test/AlgorithmsSpec.jl>>[init]
@testset "CarboKitten.Algorithms" begin
    using CarboKitten.Algorithms.RangeFinder

    # ~/~ begin <<docs/src/algorithms/enumerate_seq.md#algorithms-spec>>[init]
    a = [[:a, :b], [:c]]
    @test collect(enumerate_seq(a)) == [[(1, :a), (2, :b)], [(3, :c)]]
    # ~/~ end
end
# ~/~ end
