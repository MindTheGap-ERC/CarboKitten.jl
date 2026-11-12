# ~/~ begin <<docs/src/models/alcap.md#test/Models/ALCAPSpec.jl>>[init]
module ALCAPSpec

using Test
using CarboKitten
using CarboKitten.Models: ALCAP

@testset "Models/ALCAP/allocations" begin
    input = ALCAP.Example.INPUT
    state = ALCAP.initial_state(input)
    step! = ALCAP.step!(input)
    for _ in 1:5
        step!(state)
    end
    @test (@allocated step!(state)) < 1024
end

end
# ~/~ end
