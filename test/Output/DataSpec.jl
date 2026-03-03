# ~/~ begin <<docs/src/output/memory-writer.md#test/Output/DataSpec.jl>>[init]
module OutputDataSpec
using CarboKitten
using CarboKitten.Output.Abstract: frame_writer, add_data_set, Frame
using CarboKitten.Output.MemoryWriter: MemoryOutput
using Unitful
using Test

const DummyFacies = [
    ALCAP.Facies(
        viability_range = (0, 0),
        activation_range = (0, 0),
        maximum_growth_rate=0.0u"m/Myr",
        extinction_coefficient=0.0u"m^-1",
        saturation_intensity=0.0u"W/m^2",
        diffusion_coefficient=0.0u"m/yr"),
    ALCAP.Facies(
        viability_range = (0, 0),
        activation_range = (0, 0),
        maximum_growth_rate=0.0u"m/Myr",
        extinction_coefficient=0.0u"m^-1",
        saturation_intensity=0.0u"W/m^2",
        diffusion_coefficient=0.0u"m/yr",
        initial_sediment=3.0u"m")    
        ]

const input = ALCAP.Input(
    tag="test",
    box=Box{Periodic{2}}(grid_size=(5, 1), phys_scale=5.0u"m"),
    time=TimeProperties(
        Δt=0.0001u"Myr",
        steps=10),
    output=Dict(
        :wi1 => OutputSpec(slice=(:,:), write_interval=1),
        :wi2 => OutputSpec(slice=(:,:), write_interval=2),
        :wi3 => OutputSpec(slice=(:,:), write_interval=3),
        :wi4 => OutputSpec(slice=(:,:), write_interval=4)),
    ca_interval=1,
    initial_topography=(x, y) -> -0.0u"m",
    sea_level = t -> 0.0u"m",
    subsidence_rate=0.0u"m/Myr",
    disintegration_rate=0.0u"m/Myr",
    insolation=0.0u"W/m^2",
    sediment_buffer_size=0,
    depositional_resolution=0.0u"m",
    facies=DummyFacies)

@testset "OutputData" begin

    out = MemoryOutput(input)
    write_frame = frame_writer(input, out)

    for (k, v) in input.output
        add_data_set(out, k, v)
    end

    # create a frame of ones to be the deposition etc. each time step
    dummy_data = ones(Float64, 1, 5, 1) * u"m"
    inc = Frame(
        production = dummy_data,
        deposition = dummy_data,
        disintegration = dummy_data
    )

    write_frame(1, ALCAP.initial_frame(input))
    for t = 1:input.time.steps
        write_frame(t+1, inc)
    end

    @testset "size of output array" begin
        @test size(out.data_volumes[:wi1].deposition)[4] == 10 + 1
        @test size(out.data_volumes[:wi2].deposition)[4] == 5 + 1
        @test size(out.data_volumes[:wi3].deposition)[4] == 3 + 1
        @test size(out.data_volumes[:wi4].deposition)[4] == 2 + 1
    end

    @testset "frame written only every write_interval" begin
        @test all(out.data_volumes[:wi1].deposition[1,:,:,2:end] .≈ out.data_volumes[:wi1].write_interval*1.0u"m")
        @test all(out.data_volumes[:wi2].deposition[1,:,:,2:end] .≈ out.data_volumes[:wi2].write_interval*1.0u"m")
        @test all(out.data_volumes[:wi3].deposition[1,:,:,2:end] .≈ out.data_volumes[:wi3].write_interval*1.0u"m")
        @test all(out.data_volumes[:wi4].deposition[1,:,:,2:end] .≈ out.data_volumes[:wi4].write_interval*1.0u"m")
    end

    @testset "initial sediment frame is left empty when undefined" begin
        @test all(out.data_volumes[:wi1].deposition[1,:,:,1] .≈ 0.0u"m")
        @test all(out.data_volumes[:wi2].deposition[1,:,:,1] .≈ 0.0u"m")
        @test all(out.data_volumes[:wi3].deposition[1,:,:,1] .≈ 0.0u"m")
        @test all(out.data_volumes[:wi4].deposition[1,:,:,1] .≈ 0.0u"m")
    end

    @testset "initial sediment frame is present when defined" begin
        @test all(out.data_volumes[:wi1].deposition[2,:,:,1] .≈ 3.0u"m")
        @test all(out.data_volumes[:wi2].deposition[2,:,:,1] .≈ 3.0u"m")
        @test all(out.data_volumes[:wi3].deposition[2,:,:,1] .≈ 3.0u"m")
        @test all(out.data_volumes[:wi4].deposition[2,:,:,1] .≈ 3.0u"m")
    end

end

end
# ~/~ end
