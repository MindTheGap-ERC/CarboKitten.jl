using CarboKitten
using Test

@testset "Unit test sanity" begin
    @test 1 + 1 == 2
end

@testset "CarboKitten" begin
    using CarboKitten.SedimentStack
    include("SedimentStackSpec.jl")
    include("TransportSpec.jl")
end
