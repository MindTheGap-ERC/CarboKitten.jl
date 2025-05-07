# ~/~ begin <<docs/src/active-layer-transport.md#examples/transport/1d-advection.jl>>[init]
#| creates: docs/src/_fig/1d-advection.svg
#| collect: figures

module Script
    using CarboKitten
    using CarboKitten.Testing: transport_test_input
    using CairoMakie

    # ~/~ begin <<docs/src/active-layer-transport.md#plot-1d-evolution>>[init]
    using Printf: @sprintf
    using Unitful: ustrip

    function plot_1d_evolution!(ax::Axis, input, every=100)
    	y_idx = 1
    	(x, y) = box_axes(input.box)
    	state = ALCAP.initial_state(input)

    	plot_state() = begin
    		t = state.step * input.time.Δt
    		η = input.initial_topography.(x, y') .+ 
                state.sediment_height .-
                input.subsidence_rate * t
    		lines!(ax, x |> in_units_of(u"km"), η[:, y_idx] |> in_units_of(u"m"),
                   label=@sprintf("%.3f Myr", ustrip(t)))
    	end

    	plot_state()
    	run_model(Model{ALCAP}, input, state) do i, _
    		if mod(i, every) == 0
    			plot_state()
    		end
    	end
    end

    function plot_1d_evolution(input, every=100)
    	fig = Figure()
    	ax = Axis(fig[1, 1], xlabel="x (km)", ylabel="η (m)")

        plot_1d_evolution!(ax, input, every)

    	Legend(fig[1, 2], ax)

    	fig
    end
    # ~/~ end

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
