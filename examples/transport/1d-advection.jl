# ~/~ begin <<docs/src/active-layer-transport.md#examples/transport/1d-advection.jl>>[init]
#| creates: docs/src/_fig/1d-advection.svg
#| collect: figures

module Script
    using CarboKitten
    using CarboKitten.Testing: transport_test_input
    using CairoMakie

    include("plot-1d-evolution.jl")

    function gaussian_initial_sediment(x, y)
        exp(-(x-10u"km")^2 / (2 * (0.5u"km")^2)) * 30.0u"m"
    end

    v_const(v_max) = _ -> (Vec2(v_max, 0.0u"m/yr"), Vec2(0.0u"1/yr", 0.0u"1/yr"))

    function main()
        CairoMakie.activate!()
        input = transport_test_input(
            initial_topography = (x, y)  -> -30.0u"m",
            initial_sediment = gaussian_initial_sediment,
            disintegration_rate = 50000.0u"m/Myr",
            wave_velocity = v_const(-5u"km/Myr")
        )

        fig = plot_1d_evolution(input, 250)
        save("docs/src/_fig/1d-advection.svg", fig)
    end
end

Script.main()
# ~/~ end
