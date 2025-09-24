# ~/~ begin <<docs/src/memory-writer.md#test/OutputDataSpec.jl>>[init]
module OutputDataSpec
using CarboKitten
using CarboKitten.MemoryWriter: add_data_set
using CarboKitten.OutputData: frame_writer
using CarboKitten.Models: Frame
using Unitful
using Test

const DummyFacies = [
    ALCAP.Facies(
        viability_range = (0, 0),
        activation_range = (0, 0),
        maximum_growth_rate=0.0u"m/Myr",
        extinction_coefficient=0.0u"m^-1",
        saturation_intensity=0.0u"W/m^2",
        diffusion_coefficient=0.0u"m/yr")]

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
    inc = zeros(Frame, input)
    inc.production .= 1.0u"m"
    inc.deposition .= 1.0u"m"
    inc.disintegration .= 1.0u"m"

    for t = 1:input.time.steps
        write_frame(t, inc)
    end

    @testset "size of output array" begin
        @test size(out.data_volumes[:wi1].deposition)[4] == 10
        @test size(out.data_volumes[:wi2].deposition)[4] == 5
        @test size(out.data_volumes[:wi3].deposition)[4] == 3
        @test size(out.data_volumes[:wi4].deposition)[4] == 2
    end

    @testset "frame written only every write_interval" begin
        @test all(out.data_volumes[:wi1].deposition .≈ out.data_volumes[:wi1].write_interval*1.0u"m")
        @test all(out.data_volumes[:wi2].deposition .≈ out.data_volumes[:wi2].write_interval*1.0u"m")
        @test all(out.data_volumes[:wi3].deposition .≈ out.data_volumes[:wi3].write_interval*1.0u"m")
        @test all(out.data_volumes[:wi4].deposition .≈ out.data_volumes[:wi4].write_interval*1.0u"m")
    end

end

end
# ~/~ end