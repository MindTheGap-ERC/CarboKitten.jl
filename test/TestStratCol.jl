module StratColSpec

using Test

using Unitful

function stratigraphic_column(deposition, disintegration)
    n_times = length(deposition) - 1
    sc = zeros(typeof(1.0u"m"), n_times)

    for ts = 1:n_times
        acc = deposition[ts] - disintegration[ts]
        if acc > 0.0u"m"
            sc[ts] = acc
            continue
        end
        ts_down = ts - 1
        while acc < 0.0u"m"
            ts_down < 1 && break
            if -acc < sc[ts_down]
                sc[ts_down] -= acc
                break
            else
                acc += sc[ts_down]
                sc[ts_down] = 0.0u"m"
            end
            ts_down -= 1
        end
    end

    sc
end

@testset "Artificial_Curve" begin
    deposition = 10 .* sin.(0.1 * (1:100)) .*1.0u"m" .+ 10 .* 1.0u"m"
    disintegration = 1 .* cos.(0.1 * (1:100)).*1.0u"m" .+ 5 .* 1.0u"m"
    result = stratigraphic_column(deposition, disintegration)
    @test sum(result) ≈ sum(deposition) - sum(disintegration)
end

@testset "Strat_Col_Test" begin
    deposition = 15 * rand(20,3).*1.0u"m"
    disintegration = vcat(2 * rand(10,3), 10 * rand(10,3)).*1.0u"m"
    expected_result = deposition .- disintegration
    expected_result_first = expected_result[1:10,:]
    expected_result_sum = sum(deposition .- disintegration)
    result = stratigraphic_column(deposition, disintegration)
    @test sum(result) ≈ expected_result_sum
    @test sum(deposition) > sum(disintegration)
    @test abs(sum(expected_result_first[expected_result_first .< 0.0u"m"])) ≈ abs(sum(result) - expected_result_sum)
end

end