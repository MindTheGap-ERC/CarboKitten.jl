# ~/~ begin <<docs/src/algorithms/finite-difference-transport.md#test/Transport/AdvectionSpec.jl>>[init]
using CarboKitten
using CarboKitten: Box, box_axes
using CarboKitten.Components.TimeIntegration: time
using CarboKitten.Transport.Solvers: runge_kutta_4
using CarboKitten.Transport.Advection: transport
using CarboKitten.Testing: transport_test_input
using Unitful

@testset "CarboKitten.Transport.Advection.scale-invariance" begin

# test transport code for scale invariance
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

# ~/~ begin <<docs/src/algorithms/finite-difference-transport.md#test-wave-transport>>[init]
@testset "CarboKitten.Transport.Advection.wave-transport" begin

function gaussian_initial_sediment(x, y)
	exp(-(x-10u"km")^2 / (2 * (0.5u"km")^2)) * 30.0u"m"
end

v_const(v_max) = _ -> (Vec2(v_max, 0.0u"m/yr"), Vec2(0.0u"1/yr", 0.0u"1/yr"))

input = transport_test_input(
    initial_topography = (x, y)  -> -35.0u"m",
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

(x, y) = box_axes(input.box)

com_positions = [center_of_mass(r, x) for r in result]

@assert length(com_positions) == length(times)
observed_speeds = (com_positions[2:end] .- com_positions[1:end-1]) ./ 
    (times[2:end] .- times[1:end-1])

# Extract the expected speed from the wave transport
expected_speed = input.facies[1].wave_velocity(0.0u"m")[1][1]  
tolerance = 1e-6 * expected_speed  # Set a tolerance for the comparison
# Check if all speeds are within tolerance
for speed in observed_speeds
    @test speed ≈ expected_speed atol=tolerance
end

end
# ~/~ end
# ~/~ end
