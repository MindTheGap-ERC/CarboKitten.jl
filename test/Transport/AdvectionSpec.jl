# ~/~ begin <<docs/src/finite-difference-transport.md#test/Transport/AdvectionSpec.jl>>[init]
@testset "CarboKitten.Transport.Advection.scale-invariance" begin

using CarboKitten
using CarboKitten: Box, box_axes
using CarboKitten.Components.TimeIntegration: time
using CarboKitten.Transport.Solvers: runge_kutta_4
using CarboKitten.Transport.Advection: transport
using CarboKitten.Testing: transport_test_input
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
end

@testset "CarboKitten.Transport.Advection.advection" begin

function gaussian_initial_sediment(x, y)
	exp(-(x-10u"km")^2 / (2 * (0.5u"km")^2)) * 30.0u"m"
end

v_const(v_max) = _ -> (Vec2(v_max, 0.0u"m/yr"), Vec2(0.0u"1/yr", 0.0u"1/yr"))

input = transport_test_input(
    initial_topography = (x, y)  -> -30.0u"m",
    initial_sediment = gaussian_initial_sediment,
    disintegration_rate = 50000.0u"m/Myr",
    wave_velocity = v_const(-5u"km/Myr")
)

function center_of_mass(m, x)
    total_mass = sum(m)        # Total mass

    if total_mass ≈ 0.0u"m"
        return NaN * u"m"  # Avoid division by zero
    end

    x_com = sum(x .* sum(m, dims=2)) / total_mass  # Weighted average along x
    return x_com
end
state = ALCAP.initial_state(input)

result = []
times = []

run_model(Model{ALCAP}, input, state) do i, delta
    if mod(i-1, 250) == 0
        push!(result, deepcopy(state.sediment_height))
        push!(times, time(input, state))
    end
end
println
println("This is the result: ", result)


(x, y) = box_axes(input.box)

# Print results
# for (i, r) in enumerate(result)
#     @info "Time: $(r[2].time), Center of Mass: $(centers_of_mass[i])"
# end

fig = Figure()
ax2 = Axis(fig[1, 1], xlabel="x (km)", ylabel="η (m)")

for (i, r) in enumerate(result)
    η = input.initial_topography.(x, y') .+ r .- input.subsidence_rate * times[i]

    lines!(ax2, x |> in_units_of(u"km"), η[:, 1] |> in_units_of(u"m"),
    label=@sprintf("%.3f Myr", ustrip(times[i])))
end

com_positions = [center_of_mass(r, x) for r in result]
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

println()
println("These are the com_positions: ", com_positions)
println("and these are the times: ", times)
@assert length(com_positions) == length(times)
observed_speeds = (com_positions[2:end] .- com_positions[1:end-1]) ./ 
    (times[2:end] .- times[1:end-1])
println()
println("These are the observed speeds: ", observed_speeds)
# @show input.wave_transport             
expected_speed = input.facies[1].wave_velocity(0.0u"m")[1][1]  # Extract the expected speed from the wave transport
@show expected_speed
tolerance = 1e-6 * expected_speed  # Set a tolerance for the comparison
# Check if all speeds are within tolerance
for speed in observed_speeds
    @test speed ≈ expected_speed atol=tolerance
end
fig

end
# ~/~ end
