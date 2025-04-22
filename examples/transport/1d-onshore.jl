# ~/~ begin <<docs/src/active-layer-transport.md#examples/transport/1d-onshore.jl>>[init]
#| creates: docs/src/_fig/1d-onshore.svg
#| collect: figures

module Script
    using CarboKitten
    using CairoMakie

    # ~/~ begin <<docs/src/active-layer-transport.md#transport-test-input>>[init]
    transport_test_input(;
    	initial_topography = (x, y) -> 0.0u"m",
    	initial_sediment = (x, y) -> 0.0u"m",
    	disintegration_rate = 50.0u"m/Myr",
    	subsidence_rate = 50.0u"m/Myr",
    	diffusion_coefficient = 0.0u"m/yr",
    	wave_velocity = _ -> (Vec2(0.0, 0.0)u"m/yr", Vec2(0.0, 0.0)u"1/yr")) =

    	ALCAP.Input(
    		box = CarboKitten.Box{Coast}(grid_size=(100, 1), phys_scale=150.0u"m"),
    		time = TimeProperties(
    			Δt = 0.001u"Myr",
    			steps = 1000),
    		facies = [ALCAP.Facies(
    			initial_sediment = initial_sediment,
    			diffusion_coefficient = diffusion_coefficient,
    			wave_velocity = wave_velocity,
    			maximum_growth_rate = 0.0u"m/Myr",
    			extinction_coefficient = 0.8u"m^-1",
    			saturation_intensity = 60u"W/m^2"
    		)],
    		disintegration_rate = disintegration_rate,
    		initial_topography = initial_topography,
    		insolation = 400.0u"W/m^2",
    		sediment_buffer_size = 5,
    		depositional_resolution = 1000.0u"m",
    		transport_solver = Val{:forward_euler})
    # ~/~ end
    # ~/~ begin <<docs/src/active-layer-transport.md#plot-1d-evolution>>[init]
    using Printf: @sprintf
    using Unitful: ustrip

    function plot_1d_evolution!(ax::Axis, input, every=100)
    	y_idx = 1
    	(x, y) = box_axes(input.box)
    	state = ALCAP.initial_state(input)

    	plot_state() = begin
    		t = state.step * input.time.Δt
    		η = input.initial_topography.(x, y') .+ state.sediment_height .- input.subsidence_rate * t
    		lines!(ax, x |> in_units_of(u"km"), η[:, y_idx] |> in_units_of(u"m"), label=@sprintf("%.3f Myr", ustrip(t)))
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
