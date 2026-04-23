# ~/~ begin <<docs/src/active-layer-transport.md#test/Transport/IntertidalZoneSpec.jl>>[init]
@testset "CarboKitten.Transport.IntertidalZone" begin
    using CarboKitten
    using CarboKitten.Testing: transport_test_input

    function end_sediment_height(input)
        state = ALCAP.initial_state(input)
        run_model((_, _) -> (), Model{ALCAP}, input, state)
        return state.sediment_height
    end

    function three_peaks(x, y)
        sum(exp(-(x-μ)^2 / (2 * (0.5u"km")^2)) * 9.0u"m"
            for μ in [2.5u"km", 7.5u"km", 12.5u"km"])
    end

    function staircase(dx, dy, y0)
        return function (x, _)
            floor(x / dx) * dy + y0
        end
    end

    v_const(v_max) = _ -> (Vec2(v_max, 0.0u"m/yr"), Vec2(0.0u"1/yr", 0.0u"1/yr"))

    input1 = transport_test_input(
        initial_topography = staircase(5.0u"km", -10.0u"m", 10.0u"m"),
        initial_sediment = three_peaks,
        disintegration_rate = 50000.0u"m/Myr",
        wave_velocity = v_const(-1u"km/Myr"),
        intertidal_zone = 0u"m"
    )

    output1 = end_sediment_height(input1)[:, 1]

    input2 = transport_test_input(
        initial_topography = staircase(5.0u"km", -10.0u"m", 10.0u"m"),
        initial_sediment = three_peaks,
        disintegration_rate = 50000.0u"m/Myr",
        wave_velocity = v_const(-1u"km/Myr"),
        intertidal_zone = 10u"m"
    )

    output2 = end_sediment_height(input2)[:, 1]

    @test output1[10:30] ≈ output1[50:70] atol=0.01u"m"
    @test !isapprox(output1[50:70], output1[90:110], atol=1.0u"m")

    @test !isapprox(output2[10:30], output2[50:70], atol=1.0u"m")
    @test output2[50:70] ≈ output2[90:110] atol=0.01u"m"
end
# ~/~ end