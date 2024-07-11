using CarboKitten
using Test
using Unitful

@testset "Unit test sanity" begin
    @test 1 + 1 == 2
end

include("Unitful.jl")
include("DSLSpec.jl")

@testset "CarboKitten" begin
    include("ConfigSpec.jl")
    include("CASpec.jl")
    include("SedimentStackSpec.jl")
end
