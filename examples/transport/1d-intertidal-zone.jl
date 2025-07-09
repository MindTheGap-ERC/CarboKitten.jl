# ~/~ begin <<docs/src/active-layer-transport.md#examples/transport/1d-intertidal-zone.jl>>[init]
#| creates: docs/src/_fig/1d-intertidal.svg
#| collect: figures

module Script
    using CarboKitten
    using CarboKitten.Testing: transport_test_input
    using CairoMakie

    include("plot-1d-evolution.jl")

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

    function main()
        CairoMakie.activate!()

        input = transport_test_input(
            initial_topography = staircase(5.0u"km", -10.0u"m", 10.0u"m"),
            initial_sediment = three_peaks,
            disintegration_rate = 50000.0u"m/Myr",
            wave_velocity = v_const(-1u"km/Myr"),
            intertidal_zone = 10u"m"
        )

        fig = plot_1d_evolution(input, 500)
        save("docs/src/_fig/1d-intertidal.svg", fig)
    end
end

Script.main()
# ~/~ end
