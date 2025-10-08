# ~/~ begin <<docs/src/components/cellular-automata.md#test/Components/CellularAutomatonSpec.jl>>[init]
module CellularAutomatonSpec

using Test
using CarboKitten.Components.Common
using CarboKitten.Components: CellularAutomaton as CA

@testset "Components/CellularAutomaton" begin
    let facies = fill(CA.Facies((4, 10), (6, 10), true), 3),
        input1 = CA.Input(
            box=Box{Periodic{2}}(grid_size=(50, 50), phys_scale=1.0u"m"),
            facies=facies),
        input2 = CA.Input(
            box=Box{Periodic{2}}(grid_size=(50, 50), phys_scale=1.0u"m"),
            facies=facies,
            ca_random_seed=1)

        state1 = CA.initial_state(input1)
        state2 = CA.initial_state(input2)
        state3 = CA.initial_state(input2)

        @test CA.initial_state(input1).ca == CA.initial_state(input1).ca
        @test state1.ca != state2.ca

        step! = CA.step!(input1)  # inputs have same rules
        for i in 1:20
            step!(state1)
            step!(state2)
            step!(state3)
        end

        @test state1.ca != state2.ca
        @test state2.ca == state3.ca
    end

    @testset "inactive facies"
    let facies = [CA.Facies((4, 10), (6, 10), true),
                  CA.Facies((4, 10), (6, 10), false),
                  CA.Facies((4, 10), (6, 10), true),
                  CA.Facies((4, 10), (6, 10), false)],   
        input = CA.Input(
            box=Box{Periodic{2}}(grid_size=(10, 10), phys_scale=1.0u"m"),
            facies=facies)

        state = CA.initial_state(input)
        @test state.ca_priority == [1,3]

        # inactive facies shouldn't appear in the ca
        
        @test all(state.ca .!= 2)
        @test all(state.ca .!= 4)

    end
end

end
# ~/~ end