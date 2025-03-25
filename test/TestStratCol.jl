module StratColSpec

using Test
using CarboKitten.Export: stratigraphic_column

@testset "Strat_Col_Test" begin
    deposition = 10 * rand(10,3) .* 1.0u"m"
    disintegration = 4 * rand(10,3) .* 1.0u"m"
    expected_result = deposition .- disintegration
    result = stratigraphic_column(deposition, disintegration)
    @test result ≈ expected_result
    @test sum(deposition) > sum(disintegration)
end

end