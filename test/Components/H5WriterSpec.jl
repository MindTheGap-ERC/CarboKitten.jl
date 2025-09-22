# ~/~ begin <<docs/src/components/hdf5.md#test/Components/H5WriterSpec.jl>>[init]
module H5WriterSpec
using CarboKitten
using HDF5
using CarboKitten.Components.H5Writer: write_frame, create_ck_group
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

const filename = "testH5.h5"


@testset "Components/H5Writer" begin

    mktempdir() do path
        fpath = joinpath(path, filename)
        h5open(fpath, "w") do fid
            create_group(fid, "input")
            ALCAP.write_header(fid, input)

            for (k, v) in input.output
                create_ck_group(fid, input, k, v)
            end

            # create a frame of ones to be the deposition etc. each time step
            inc = zeros(Frame, input)
            inc.production .= 1.0u"m"
            inc.deposition .= 1.0u"m"
            inc.disintegration .= 1.0u"m"

            for t = 1:input.time.steps
                write_frame(fid, input, t, inc)
            end
        end

        @testset "size of output array" begin
            h5open(fpath, "r") do f
                @test size(f["wi1"]["deposition"][])[4] == 10
                @test size(f["wi2"]["deposition"][])[4] == 5
                @test size(f["wi3"]["deposition"][])[4] == 3
                @test size(f["wi4"]["deposition"][])[4] == 2
            end
        end

        @testset "frame written only every write_interval" begin
            h5open(fpath, "r") do f
                @test all(f["wi1"]["deposition"][] .≈ attrs(f["wi1"])["write_interval"])
                @test all(f["wi2"]["deposition"][] .≈ attrs(f["wi2"])["write_interval"])
                @test all(f["wi3"]["deposition"][] .≈ attrs(f["wi3"])["write_interval"])
                @test all(f["wi4"]["deposition"][] .≈ attrs(f["wi4"])["write_interval"])
            end
        end
    end
end

end
# ~/~ end
