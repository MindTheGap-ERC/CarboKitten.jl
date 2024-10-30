# ~/~ begin <<docs/src/utility.md#test/UtilitySpec.jl>>[init]
@testset "CarboKitten.Utility" begin
    using CarboKitten.Utility

    # ~/~ begin <<docs/src/utility.md#utility-spec>>[init]
    a = [[:a, :b], [:c]]
    @test collect(enumerate_seq(a)) == [[(1, :a), (2, :b)], [(3, :c)]]
    # ~/~ end
end
# ~/~ end
