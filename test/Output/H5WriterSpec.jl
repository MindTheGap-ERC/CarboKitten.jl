# ~/~ begin <<docs/src/output/h5writer.md#test/Output/H5WriterSpec.jl>>[init]
module H5WriterSpec

using CarboKitten
using HDF5
using CarboKitten.Output.Abstract: frame_writer, add_data_set, Frame
using CarboKitten.Output.H5Writer: H5Output
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
        initial_sediment=3.0u"m"),
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

const filename = "testH5.h5"


@testset "Components/H5Writer" begin
    mktempdir() do path
        fpath = joinpath(path, filename)
        output = H5Output(input, fpath)
        write_frame = frame_writer(input, output)

        ALCAP.write_header(input, output)

        for (k, v) in input.output
            add_data_set(output, k, v)
        end

        # create a frame of ones to be the deposition etc. each time step
        dummy_data = ones(Float64, 2, 5, 1) * u"m"
        inc = Frame(
            production = dummy_data,
            deposition = dummy_data,
            disintegration = dummy_data
        )

        write_frame(1, ALCAP.initial_frame(input))
        for t = 1:input.time.steps
            write_frame(t+1, inc)
        end

        close(output.fid)

        @testset "size of output array" begin
            h5open(fpath, "r") do f
                @test size(f["wi1"]["deposition"][])[4] == 10 + 1
                @test size(f["wi2"]["deposition"][])[4] == 5 + 1
                @test size(f["wi3"]["deposition"][])[4] == 3 + 1
                @test size(f["wi4"]["deposition"][])[4] == 2 + 1
            end
        end

        @testset "frame written only every write_interval" begin
            h5open(fpath, "r") do f
                @test all(f["wi1"]["deposition"][][1,:,:,2:end] .≈ attrs(f["wi1"])["write_interval"])
                @test all(f["wi2"]["deposition"][][1,:,:,2:end] .≈ attrs(f["wi2"])["write_interval"])
                @test all(f["wi3"]["deposition"][][1,:,:,2:end] .≈ attrs(f["wi3"])["write_interval"])
                @test all(f["wi4"]["deposition"][][1,:,:,2:end] .≈ attrs(f["wi4"])["write_interval"])
            end
        end

        @testset "initial sediment frame is left empty when undefined" begin
            h5open(fpath, "r") do f
                @test all(f["wi1"]["deposition"][][1,:,:,1] .≈ 0.0)
                @test all(f["wi2"]["deposition"][][1,:,:,1] .≈ 0.0)
                @test all(f["wi3"]["deposition"][][1,:,:,1] .≈ 0.0)
                @test all(f["wi4"]["deposition"][][1,:,:,1] .≈ 0.0)
            end
        end

        @testset "initial sediment frame is present when defined" begin
            h5open(fpath, "r") do f
                @test all(f["wi1"]["deposition"][][2,:,:,1] .≈ 3.0)
                @test all(f["wi2"]["deposition"][][2,:,:,1] .≈ 3.0)
                @test all(f["wi3"]["deposition"][][2,:,:,1] .≈ 3.0)
                @test all(f["wi4"]["deposition"][][2,:,:,1] .≈ 3.0)
            end
        end

    end
end

end
# ~/~ end
