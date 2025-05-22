module StratColSpec

using Test

using Unitful

function stratigraphic_column(deposition, disintegration)
    @assert size(deposition) == size(disintegration)

    n_facies = size(deposition)[2]
    n_times = size(deposition)[1]

    sc = zeros(typeof(1.0u"m"), n_times, n_facies)

    for ts = 1:n_times
        acc = deposition[ts, :] .- disintegration[ts, :]
        @assert acc <= deposition[ts, :]
        if sum(acc) > 0.0u"m"
            sc[ts, :] .= acc
            @assert sum(sc[ts, :]) ≈ sum(acc)
            @assert sum(sc[ts, :]) > 0.0u"m"
            continue
        end

        ts_down = ts - 1
        @assert ts_down >= 1
        while sum(acc) < 0.0u"m"
            if ts_down < 0
                @warn "stratigraph column overshoot: $(acc)"
                break
            end
            ts_down < 1 && break
            if -sum(acc) < sum(sc[ts_down, :])
                previous_sc = sc[ts_down, :]
                sc[ts_down, :] .-= acc
                @assert sum(sc[ts_down, :]) >= 0.0u"m"
                @assert sum(sc[ts_down, :]) + sum(acc) ≈ sum(previous_sc)
                if any(sc[ts_down, :] .< 0.0u"m")
                    @warn "negative value in stratigraphic column: $(sc[ts_down,:])"
                end
                break
            end
            
            acc .+= sc[ts_down, :]
            @assert acc .<= 0.0u"m"
            if any(acc .> 0.0u"m")
                @warn "round-off error in stratigraphic column: $(acc)"
            end
            sc[ts_down] .= 0.0u"m"
            ts_down -= 1
        end
    end

    return sc
end

@testset "Artificial_Curve" begin
    deposition = 10 .* sin.(0.1 * (1:100)) .*1.0u"m" .+ 10 .* 1.0u"m"
    disintegration = 1 .* sin.(0.1 * (1:100)).*1.0u"m" .+ 5 .* 1.0u"m"
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