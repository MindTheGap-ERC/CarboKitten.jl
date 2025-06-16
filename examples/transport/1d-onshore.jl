# ~/~ begin <<docs/src/active-layer-transport.md#examples/transport/1d-onshore.jl>>[init]
#| creates: docs/src/_fig/1d-onshore.svg
#| collect: figures

module Script
    using CarboKitten
    using CarboKitten.Testing: transport_test_input
    using CairoMakie

    include("plot-1d-evolution.jl")

    v_prof(v_max, max_depth, w) = 
        let k = sqrt(0.5) / max_depth,
            A = 3.331 * v_max,
            α = tanh(k * w),
            β = exp(-k * w)
            (A * α * β, -A * k * β * (1 - α - α^2))
        end

    v_const(v_max) = _ -> (Vec2(v_max, 0.0u"m/yr"), Vec2(0.0u"1/yr", 0.0u"1/yr"))

    v_prof_par(v_max, max_depth) = w -> let (v, s) = v_prof(v_max, max_depth, w)
        (Vec2(v, 0.0u"m/yr"), Vec2(s, 0.0u"1/yr"))
    end

    function main()
        CairoMakie.activate!()

        fig = Figure(size=(1000, 500))
        ax1 = Axis(fig[1, 1], xlabel="x (km)", ylabel="η (m)")
        input1 = transport_test_input(
            initial_topography = (x, y) -> -x / 375.0 - 10u"m",
            initial_sediment = 10.0u"m",
            disintegration_rate = 50.0u"m/Myr",
            diffusion_coefficient = 5.0u"m/yr",
            wave_velocity = v_const(-0.5u"m/yr")
        )
        plot_1d_evolution!(ax1, input1, 250)

        ax2 = Axis(fig[1, 2], xlabel="x (km)", ylabel="η (m)")
        input2 = transport_test_input(
            initial_topography = (x, y) -> -x / 375.0 - 10u"m",
            initial_sediment = 10.0u"m",
            disintegration_rate = 50.0u"m/Myr",
            diffusion_coefficient = 5.0u"m/yr",
            wave_velocity = v_prof_par(-0.5u"m/yr", 20u"m")
        )
        plot_1d_evolution!(ax2, input2, 250)

        Legend(fig[1, 3], ax2)
        save("docs/src/_fig/1d-onshore.svg", fig)
    end
end

Script.main()
# ~/~ end
