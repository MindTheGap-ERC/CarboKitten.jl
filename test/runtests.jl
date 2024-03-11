using CarboKitten
using Test
using Unitful

@testset "Unit test sanity" begin
    @test 1 + 1 == 2
end

include("Unitful.jl")

@testset "CarboKitten" begin
    using CarboKitten.SedimentStack
    include("SedimentStackSpec.jl")
    include("TransportSpec.jl")
end
