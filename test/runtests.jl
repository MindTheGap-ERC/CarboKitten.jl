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
    include("DenudationSpec.jl")
    include("SedimentStackSpec.jl")

    include("Components/TimeIntegrationSpec.jl")
    include("Components/FaciesBaseSpec.jl")
    include("Components/TimeIntegrationSpec.jl")
    include("Components/BoxesSpec.jl")
    include("Components/ProductionSpec.jl")

    include("ExportSpec.jl")
end
