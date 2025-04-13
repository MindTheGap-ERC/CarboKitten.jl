# Initial Sediment

```component-dag
CarboKitten.Components.InitialSediment
```

``` {.julia file=src/Components/InitialSediment.jl}
@compose module InitialSediment
    @mixin FaciesBase, WaterDepth, SedimentBuffer
    using Unitful: dimension, NoUnits

    <<initial-sediment>>
end
```

Mostly for debugging and testing purposes it can be convenient to specify an initial layer of sediment. This component adds the `initial_sediment` parameter to the facies definition.

``` {.julia #initial-sediment}
@kwdef struct Facies <: AbstractFacies
    initial_sediment
end
```

If there are multiple facies that define an initial sediment, those are all added to the sediment buffer in one go. The sediment buffer will be filled with a uniform ratio of facies.

```  {.julia #initial-sediment}
function initial_sediment(box::Box, facies::AbstractFacies)
    s = facies.initial_sediment
    if s isa AbstractMatrix
        @assert size(s) == box.grid_size
        @assert dimension(eltype(s)) == dimension(u"m")
        return s
    end
    if s isa Quantity
        @assert dimension(s) == dimension(u"m")
        return fill(s, box.grid_size...)
    end

    # s should be callable
    x, y = box_axes(box)
    return s(x, y')
end

function push_initial_sediment!(input::AbstractInput, state::AbstractState)
    s = stack(initial_sediment(input.box, f) for f in input.facies; dims=1)
    push_sediment!(state.sediment_buffer, s ./ input.depositional_resolution .|> NoUnits)
    state.sediment_height .+= sum(s; dims=1)[1,:,:]
end
```

## Tests

``` {.julia file=test/Components/InitialSedimentSpec.jl}
@testset "Components/InitialSediment" begin

using Unitful
using CarboKitten.Components.Common
using CarboKitten.Components: InitialSediment as IS

function initial_state(input::AbstractInput)
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, IS.n_facies(input), input.box.grid_size...)
    state = IS.State(
        step = 0,
        sediment_height = zeros(IS.Sediment, input.box.grid_size...),
        sediment_buffer = sediment_buffer)
    IS.push_initial_sediment!(input, state)
    return state
end

let box = Box{Periodic}(grid_size=(50, 50), phys_scale=100.0u"m"),
    input = IS.Input(
        time = TimeProperties(steps=0, Δt=1.0u"yr"),
        box = box,
        facies = [ IS.Facies(initial_sediment=30u"m") ] ),
    state = initial_state(input)

    @test all(isapprox.(state.sediment_height, 30.0u"m"))
    @test all(isapprox.(state.sediment_buffer, 1.0))
end

end
```
