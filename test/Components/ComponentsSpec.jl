# ~/~ begin <<docs/src/components/components.md#test/Components/ComponentsSpec.jl>>[init]
@testset "Components" begin
    using CarboKitten.Components: Components, list_components

    # every exported name exists
    @test all(n -> isdefined(Components, n), names(Components))

    components = list_components()
    @test all(c -> c isa Module, components)
    @test Components.SedimentBuffer in components
    @test Components.CircularSedimentBuffer in components
    @test Components.Common ∉ components
    # every component is exported
    @test all(c -> nameof(c) in names(Components), components)
end
# ~/~ end
