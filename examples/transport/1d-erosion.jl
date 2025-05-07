# ~/~ begin <<docs/src/active-layer-transport.md#examples/transport/1d-erosion.jl>>[init]
#| creates: docs/src/_fig/1d-erosion.svg
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

	function initial_sediment(x, y)
	  if x < 5.0u"km"
	    return 30.0u"m"
	  end

	  if x > 10.0u"km" && x < 11.0u"km"
	    return 20.0u"m"
	  end

	  return 5.0u"m"
	end

    function main()
        CairoMakie.activate!()
        input = transport_test_input(
            initial_topography = (x, y) -> -30.0u"m",
            initial_sediment = initial_sediment,
            diffusion_coefficient = 10.0u"m/yr")

        fig = plot_1d_evolution(input, 250)
        save("docs/src/_fig/1d-erosion.svg", fig)
    end
end

Script.main()
# ~/~ end
