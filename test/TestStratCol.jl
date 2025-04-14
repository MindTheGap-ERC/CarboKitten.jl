module StratColSpec

using Test
using CarboKitten.Export: stratigraphic_column
using Unitful
@testset "Strat_Col_Test" begin
    deposition = 15 * rand(20,3) .* 1.0u"m"
    disintegration = vcat(2 * rand(10,3) .* 1.0u"m", 10 * rand(10,3) .* 1.0u"m")
    expected_result = deposition .- disintegration
    expected_result_first = expected_result[1:10,:]
    expected_result_sum = sum(deposition .- disintegration)
    result = stratigraphic_column(deposition, disintegration)
    @test sum(result) ≈ expected_result_sum
    @test sum(deposition) > sum(disintegration)
    @test abs(sum(expected_result_first[expected_result_first .< 0.0u"m"])) ≈ abs(sum(result) - expected_result_sum)
end

end