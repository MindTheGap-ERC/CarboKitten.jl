# ~/~ begin <<docs/src/carbocat-ca.md#test/CASpec.jl>>[init]
@testset "CA" begin
    using CarboKitten.BoundaryTrait: Periodic
    using CarboKitten.Config: Box
    using CarboKitten.Burgess2013.CA: step_ca, run_ca
    using Unitful

    using CarboKitten.Components: CellularAutomaton as CA

    MODEL1 = [
        CA.Facies(viability_range=(4, 10), activation_range=(6, 10)),
        CA.Facies(viability_range=(4, 10), activation_range=(6, 10)),
        CA.Facies(viability_range=(4, 10), activation_range=(6, 10))]

    n_facies = length(MODEL1)
    box = Box{Periodic{2}}(grid_size=(50, 50), phys_scale=100.0u"m")
    ca_init = rand(0:n_facies, box.grid_size...)
    ca_channel = run_ca(Periodic{2}, MODEL1, copy(ca_init), n_facies)
    item1, _ = Iterators.peel(Iterators.drop(ca_channel, 20))

    ca_step = step_ca(box, MODEL1)
    state = CA.State(copy(ca_init), 1:n_facies)
    for _ in 1:20
        ca_step(state)
    end

    @test item1 == state.ca
end
# ~/~ end
