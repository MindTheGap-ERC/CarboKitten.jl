# ~/~ begin <<docs/src/components/sediment_buffer.md#test/Components/SedimentBufferSpec.jl>>[init]
module SedimentBufferSpec
using Test
using Unitful
using CarboKitten
using CarboKitten.Components.Common: Amount, record_type
using CarboKitten.SedimentBuffers: initial_sediment_buffer, peek_sediment
import CarboKitten.Components.SedimentBuffer as SB
import CarboKitten.Components: CopyingSedimentBuffer, CircularSedimentBuffer

@testset "Components/SedimentBuffer" begin
    @testset "$(impl)" for impl in (CopyingSedimentBuffer, CircularSedimentBuffer)
        input = record_type(:Input, SB, impl)(
            box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
            time = TimeProperties(Δt=1.0u"yr", steps=10),
            facies = [record_type(:Facies, SB, impl)()],
            sediment_buffer_size = 10,
            depositional_resolution = 1.0u"m")
        state = record_type(:State, SB, impl)(
            step = 0,
            bathymetry = zeros(Amount, input.box.grid_size...),
            sediment_buffer = initial_sediment_buffer(input),
            sediment_thickness = zeros(Amount, input.box.grid_size...))

        @test size(state.sediment_thickness) == (10, 1)
        @test state.sediment_buffer isa fieldtype(impl.State, :sediment_buffer)

        push! = SB.push_sediment(input)
        pop! = SB.pop_sediment(input)

        push!(state, reshape((1:10) .* 0.5u"m" |> collect, (1, 10, 1)))
        @test peek_sediment(state, 10.0)[1, :, 1] ≈ (1:10) .* 0.5
        @test peek_sediment(state, CartesianIndex(3, 1), 1.0) ≈ [1.0]
        # no initial topography, no subsidence
        @test state.sediment_thickness ≈ state.bathymetry

        buffer = zeros(Amount, 1, 10, 1)
        pop!(state, state.sediment_thickness ./ 2, buffer)
        @test reshape(state.sediment_thickness, (1, 10, 1)) ≈ buffer
        @test peek_sediment(state, 10.0)[1, :, 1] ≈ (1:10) .* 0.25
        # no initial topography, no subsidence
        @test state.sediment_thickness ≈ state.bathymetry
    end

    @testset "copying layout" begin
        input = record_type(:Input, SB)(
            box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
            time = TimeProperties(Δt=1.0u"yr", steps=10),
            facies = [record_type(:Facies, SB)()],
            sediment_buffer_size = 10,
            depositional_resolution = 1.0u"m")
        state = record_type(:State, SB)(
            step = 0,
            bathymetry = zeros(Amount, input.box.grid_size...),
            sediment_buffer = initial_sediment_buffer(input),
            sediment_thickness = zeros(Amount, input.box.grid_size...))

        @test size(state.sediment_buffer) == (10, 1, 10, 1)

        SB.push_sediment(input)(state, reshape((1:10) .* 0.5u"m" |> collect, (1, 10, 1)))
        @test state.sediment_buffer[1, 1, :, 1] == repeat([0.5, 0.0], 5)

        SB.pop_sediment(input)(state, state.sediment_thickness ./ 2, zeros(Amount, 1, 10, 1))
        @test state.sediment_buffer[1, 1, :, 1] .% 1.0 ≈ repeat([0.25, 0.5, 0.75, 0.0], 3)[1:10]
    end
end
end
# ~/~ end
