# ~/~ begin <<docs/src/finite-difference-transport.md#test/Transport/AdvectionSpec.jl>>[init]
@testset "CarboKitten.Transport.Advection.scale-invariance" begin

using CarboKitten: Box, box_axes
using CarboKitten.Transport.Solvers: runge_kutta_4
using CarboKitten.Transport.Advection: transport
using Unitful
using Printf

Amount = typeof(1.0u"m")

let box = Box{Periodic{2}}(grid_size=(32, 32), phys_scale=1.0u"m")
    solver = runge_kutta_4(Float64, box)
    wave_velocity = _ -> ((0.5u"m/s", 0.0u"m/s"), (0.0u"1/s", 0.0u"1/s"))
    diffusivity = 5.0u"m/s"
    w = randn(box.grid_size...) * u"m"
    C1 = randn(box.grid_size...)
    C2 = C1 .* 10.0
    dt = 1.0u"s"
    df(C, _) = transport(box, diffusivity, wave_velocity, C, w)

    for i = 1:10
        solver(df, C1, 0.0u"s", dt)
        solver(df, C2, 0.0u"s", dt)
    end

    @test isapprox(C1 .* 10.0, C2)
end

mutable struct State
    time::typeof(1.0u"Myr")
    sediment::Matrix{typeof(1.0u"m")}
end

function initial_state(input)
    x, y = box_axes(input.box)
    State(0.0u"Myr", input.initial_sediment.(x, y'))
end

struct Frame
    t::typeof(1.0u"Myr")
    δ::Matrix{Amount}
end

function propagator(input)
    x, y = box_axes(input.box)
    μ0 = input.initial_topography.(x, y')
    box = input.box
    Δt = input.Δt
    disintegration_rate = input.disintegration_rate
    production = input.production
    d = input.diffusivity
	v = input.wave_transport
	σ = input.subsidence_rate
	max_amount = disintegration_rate * Δt

	@info "disintegration per time step: $max_amount"
	
	solver = if input.transport_solver === nothing
		runge_kutta_4(Amount, box)
	else
		input.transport_solver
	end
	
    function active_layer(state)    
        amount = min.(max_amount, state.sediment)
        state.sediment .-= amount

        production.(x, y') * Δt .+ amount
    end

    function (state)
        p = active_layer(state)
		# @info "production between $(minimum(p)) - $(maximum(p))"
		water_depth = (σ * state.time) .- (μ0 .+ state.sediment)

		for j in 1:input.transport_substeps
			solver(
				(C, _) -> transport(box, d, v, p, water_depth),
				p, state.time, Δt / input.transport_substeps)
        end

        return Frame(state.time, ifelse.(p .< 0.0u"m", 0.0u"m", p))
    end
end

function run_model(input)
    state = initial_state(input)
    prop = propagator(input)

    Channel{State}() do ch
        while state.time < input.t_end
            Δ = prop(state)
            state.sediment .+= Δ.δ
            state.time += input.Δt
            put!(ch, state)
        end
    end
end

function gaussian_initial_sediment(x, y)
	exp(-(x-10u"km")^2 / (2 * (0.5u"km")^2)) * 30.0u"m"
end

v_const(v_max) = _ -> ((v_max, 0.0u"m/yr"), (0.0u"1/yr", 0.0u"1/yr"))

@kwdef struct Input
    box
    Δt::typeof(1.0u"Myr")
    t_end::typeof(1.0u"Myr")
    initial_topography   # function (x::u"m", y::u"m") -> u"m"
    initial_sediment    # function (x::u"m", y::u"m") -> u"m"
    production          # function (x::u"m", y::u"m") -> u"m/s"
    disintegration_rate::typeof(1.0u"m/Myr")
    subsidence_rate::typeof(1.0u"m/Myr")
    diffusivity::typeof(1.0u"m/yr")
	wave_transport
	transport_substeps::Int = 1
	transport_solver = nothing
end

# Peak should move with velocity of the onshore transport velocity.
INPUT = Input(
	box                   = Box{Coast}(grid_size=(100, 1), phys_scale=150.0u"m"),
    Δt                    = 0.01u"Myr",
    t_end                 = 1.0u"Myr",

    initial_topography     = (x, y) -> -30.0u"m",
    initial_sediment      = gaussian_initial_sediment,
    production            = (x, y) -> 0.0u"m/Myr",

    disintegration_rate   = 50000.0u"m/Myr",
    subsidence_rate       = 50.0u"m/Myr",
    diffusivity           = 0.0u"km/Myr",
	wave_transport        = v_const(-0.005u"m/yr"),
	transport_substeps    = 10
		# w -> let (v, s) = v_prof(-5u"m/yr", 20.0u"m", w)
		# 	((v, 0.0u"m/yr"), (s, 0.0u"1/yr"))
		# end
)

function center_of_mass(state, x)
    M = state.sediment         # Sediment distribution (assumed to be a 2D matrix)
    total_mass = sum(M)        # Total mass

    if total_mass ≈ 0.0u"m"
        return NaN * u"m"  # Avoid division by zero
    end

    x_com = sum(x .* sum(M, dims=2)) / total_mass  # Weighted average along x
    return x_com
end


function plot_erosion(input, every=40)
	y_idx = 1
	result = Iterators.map(deepcopy,
		Iterators.filter(x -> mod(x[1]-1, every) == 0, 
		enumerate(run_model(input)))) |> collect

	(x, y) = box_axes(input.box)
    
    # Print results
    # for (i, r) in enumerate(result)
    #     @info "Time: $(r[2].time), Center of Mass: $(centers_of_mass[i])"
    # end

	fig = Figure()
	ax2 = Axis(fig[1, 1], xlabel="x (km)", ylabel="η (m)")

	for (i, r) in enumerate(result)
		η = input.initial_topography.(x, y') .+ r[2].sediment .- 
        input.subsidence_rate * r[2].time

		lines!(ax2, x |> in_units_of(u"km"), η[:, y_idx] |> in_units_of(u"m"),
        label=@sprintf("%.3f Myr", ustrip(r[2].time)))
	end

    com_positions = [center_of_mass(r[2], x) for r in result]
    # @info "Time: $(r[2].time), Center of Mass: $(com)"
    com_positions_km = com_positions |> in_units_of(u"km")  # Convert to km
    com_positions_km_num = ustrip.(com_positions_km)  # Strip units for plotting

    # Plot centers of mass as dots
    y_values_m = fill(ustrip(u"m", -10u"m"), length(com_positions_km_num))  # Convert -10m to numeric
    # Plot centers of mass as dots
    scatter!(ax2, com_positions_km_num, y_values_m,
                 markersize=10, color=:red, label="Center of Mass")

    lines!(ax2, com_positions_km_num, y_values_m, color=:red, linestyle=:dash)

	Legend(fig[1, 2], ax2)

	display(fig)

    observed_speeds = [(com_positions[i+1] - com_positions[i]) / (result[i+1][2].time - result[i][2].time) 
                    for i in 1:length(com_positions)-1]
    @show input.wave_transport             
    expected_speed = input.wave_transport(0.0u"m")[1][1]  # Extract the expected speed from the wave transport
    @show expected_speed
    tolerance = 1e-6 * expected_speed  # Set a tolerance for the comparison
    # Check if all speeds are within tolerance
    for speed in observed_speeds
        @test speed ≈ expected_speed atol=tolerance
    end
end

plot_erosion(INPUT)

end
# ~/~ end
