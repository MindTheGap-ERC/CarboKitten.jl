module LiveView
    using CarboKitten.CaProd
    using CarboKitten.BoundaryTrait
    using CarboKitten.Config: Box, TimeProperties
    using CarboKitten.Burgess2013: production_rate, Facies
    using CarboKitten.Visualization: plot_facies_production
    using CarboKitten.Utility

    using Unitful
    using GLMakie
    using Observables
    using ProgressBars

	const PERIOD = 200.0u"kyr"
	const AMPLITUDE = 4.0u"m"
	const FACIES = [
	    Facies(viability_range = (4, 10),
	           activation_range = (6, 10),
	           maximum_growth_rate = 500u"m/Myr",
	           extinction_coefficient = 0.8u"m^-1",
	           saturation_intensity = 60u"W/m^2"),
	
	    Facies(viability_range = (4, 10),
	           activation_range = (6, 10),
	           maximum_growth_rate = 400u"m/Myr",
	           extinction_coefficient = 0.1u"m^-1",
	           saturation_intensity = 60u"W/m^2"),
	
	    Facies(viability_range = (4, 10),
	           activation_range = (6, 10),
	           maximum_growth_rate = 100u"m/Myr",
	           extinction_coefficient = 0.005u"m^-1",
	           saturation_intensity = 60u"W/m^2")
	]

	const INPUT = CaProd.Input(
	  box = Box{Shelf}(
	    grid_size = (100, 50),
	    phys_scale = 1.0u"km"
	  ),
	  time = TimeProperties(
	    Δt = 0.001u"Myr",
	    steps = 1000,
	    write_interval = 1
	  ),
	  sea_level = t -> AMPLITUDE * sin(2π * t / PERIOD), 
	  subsidence_rate=50.0u"m/Myr",
	  initial_depth=x -> x / 2000.0,
	  facies=FACIES,
	  insolation=400.0u"W/m^2"
	)

    function main(input)
        GLMakie.activate!()
        fig = Figure()
        ax = Axis3(fig[1:2, 1]; limits=(nothing, nothing, (-100, 50)))
        plot_facies_production(input; loc=fig[1,2])

        n_timesteps = input.time.steps ÷ input.time.write_interval
        height = Observable{Matrix{Float64}}(zeros(Float64, input.box.grid_size...))
        section = Observable{Matrix{Float64}}(zeros(Float64, input.box.grid_size[1], n_timesteps))
        ax_hm = Axis(fig[2,2])
        heatmap!(ax_hm, section)
        x_axis = (0:(input.box.grid_size[1]-1)) .* input.box.phys_scale
        y_axis = (0:(input.box.grid_size[2]-1)) .* input.box.phys_scale
        surface!(ax, x_axis |> in_units_of(u"m"), y_axis |> in_units_of(u"m"), height)

        function run_model(input)
            Channel{CaProd.Frame}() do ch
                s = CaProd.initial_state(input)
                p = CaProd.propagator(input)
                u = CaProd.updater(input)

                while true
                    Δ = p(s)
                    put!(ch, Δ)
                    u(s, Δ)
                    height[] = .-s.height / u"m"
                end
            end
        end

        # for f in ProgressBar(Iterators.take(run_model(INPUT), 1000), total=1000)
        
        @async begin
            for (i, f) in enumerate(Iterators.take(run_model(INPUT), 1000))
                section[][:, i] = sum(f.production[:,25,:];dims=2)
            end
        end

        return fig
    end
end  # module

LiveView.main(LiveView.INPUT)
