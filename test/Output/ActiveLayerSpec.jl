# ~/~ begin <<docs/src/active-layer-transport.md#test/Output/ActiveLayerSpec.jl>>[init]
module ActiveLayerSpec

using CarboKitten
using CarboKitten.Output.Abstract: state_writer, add_data_set, Frame
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
        diffusion_coefficient=0.0u"m/yr")]

const input = ALCAP.Input(
    tag="test",
    box=Box{Periodic{2}}(grid_size=(5, 1), phys_scale=5.0u"m"),
    time=TimeProperties(
        Δt=0.0001u"Myr",
        steps=10),
    output=Dict(
        :wi1 => OutputSpec(slice=(:,:), write_interval=1)),
    ca_interval=1,
    initial_topography=(x, y) -> -0.0u"m",
    sea_level = t -> 0.0u"m",
    subsidence_rate=0.0u"m/Myr",
    disintegration_rate=0.0u"m/Myr",
    insolation=0.0u"W/m^2",
    sediment_buffer_size=0,
    depositional_resolution=0.0u"m",
    facies=DummyFacies,
    save_active_layer=true)

@testset "ActiveLayer" begin
    out = MemoryOutput(input)
    write_state = state_writer(input, out)

    for (k, v) in input.output
        add_data_set(out, k, v)
    end

    @test out.data_volumes[:wi1].active_layer != nothing
end

end
# ~/~ end
