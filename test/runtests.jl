using CarboKitten
using Test
using Unitful

@testset "Unit test sanity" begin
    @test 1 + 1 == 2
end

include("Unitful.jl")

@testset "CarboKitten" begin
    include("ConfigSpec.jl")
    include("CASpec.jl")
    include("Denudationtest.jl")
    include("SedimentStackSpec.jl")


end
