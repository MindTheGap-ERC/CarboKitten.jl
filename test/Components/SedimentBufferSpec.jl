# ~/~ begin <<docs/src/components/sediment_buffer.md#test/Components/SedimentBufferSpec.jl>>[init]
using CarboKitten
using CarboKitten.Components.Common: Amount
import CarboKitten.Components.SedimentBuffer as SB

@testset "Components/SedimentBuffer" begin
    input = SB.Input(
        box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
        time = TimeProperties(Δt=1.0u"yr", steps=10),
        facies = [SB.Facies()],
        sediment_buffer_size = 10,
        depositional_resolution = 1.0u"m")
    state = SB._initial_state(input)

    @test size(state.sediment_thickness) == (10, 1)
    @test size(state.sediment_buffer) == (10, 1, 10, 1)

    push! = SB.push_sediment(input)
    pop! = SB.pop_sediment(input)

    push!(state, reshape((1:10) .* 0.5u"m" |> collect, (1, 10, 1)))
    @test state.sediment_buffer[1, 1, :, 1] == repeat([0.5, 0.0], 5)
    # no initial topography, no subsidence
    @test state.sediment_thickness ≈ state.bathymetry

    buffer = zeros(Amount, 1, 10, 1)
    pop!(state, state.sediment_thickness ./ 2, buffer)
    @test reshape(state.sediment_thickness, (1, 10, 1)) ≈ buffer
    @test state.sediment_buffer[1, 1, :, 1] .% 1.0 ≈ repeat([0.25, 0.5, 0.75, 0.0], 3)[1:10]
    # no initial topography, no subsidence
    @test state.sediment_thickness ≈ state.bathymetry
end
# ~/~ end
