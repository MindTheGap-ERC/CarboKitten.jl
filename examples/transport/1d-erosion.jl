# ~/~ begin <<docs/src/active-layer-transport.md#examples/transport/1d-erosion.jl>>[init]
#| creates: docs/src/_fig/1d-erosion.svg
#| collect: figures

module Script
    using CarboKitten
    using CarboKitten.Testing: transport_test_input
    using CairoMakie

    include("plot-1d-evolution.jl")

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
