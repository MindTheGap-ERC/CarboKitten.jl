# ~/~ begin <<docs/src/active-layer-transport.md#test/Components/ActiveLayerSpec.jl>>[init]
module ActiveLayerSpec

using Test
using CarboKitten.Components: ActiveLayer as AL
using CarboKitten.Components.Common

@testset "Components/ActiveLayer" begin
    @testset "Disintegration transfer" begin
        let facies = fill(AL.Facies(), 4),
            input = AL.Input(
                box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
                time = TimeProperties(Δt=1.0u"kyr", steps=10),
                facies=facies,
                disintegration_transfer = f -> stack((0.0.*f[1,:,:], 0.5.*f[2,:,:], 
                                          f[1,:,:].+f[3,:,:], f[4,:,:].+0.5.*f[2,:,:]),dims=1),
            )

            state = AL.State(
                step = 0,
                sediment_height = zeros(Height, input.box.grid_size...),
                sediment_buffer = zeros(Float64, input.sediment_buffer_size, AL.n_facies(input), input.box.grid_size...),
                active_layer=zeros(Amount, AL.n_facies(input), input.box.grid_size...))


            d = ones(Amount, AL.n_facies(input), input.box.grid_size...)

            dtf = input.disintegration_transfer
            transferred_sed = dtf(d)
            state.active_layer .+= transferred_sed 

            @test all(state.active_layer[1,:] .≈ 0.0u"m")
            @test all(state.active_layer[2,:] .≈ 0.5u"m")
            @test all(state.active_layer[3,:] .≈ 2.0u"m")
            @test all(state.active_layer[4,:] .≈ 1.5u"m")

            @test sum(state.active_layer[:,1,1]) ≈ sum(d[:,1,1])

            @test all(size(transferred_sed) .== size(d))

        end
    end
end
end
# ~/~ end
