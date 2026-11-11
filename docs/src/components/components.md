# Model Components

Each model in CarboKitten is composed out of elementary parts. These form a hierarchy of modules that inherit from each other.

<!-- FIXME: auto-generate this graph

Mermaid is broken...

-->

Not all of these are equally important. The `H5Writer` implements the core simulation loop and writing a simulation to HDF5 files. For a `Model` to work with `H5Writer` it needs to have `initial_state`, `step!`, and `write_header` functions implemented. Each component can have its own `write_header(::Input, input, output)` implementation, dispatching on its own `Input` part. The generic `write_header(input, output)` calls it for every part of the input.

The initial state is composed in the same way. Each component that has a `State` part implements `initial_state(::Input, input)`, returning its part of the state. The generic `initial_state(State, input)` collects these parts using `map_components` from ModuleMixins. When a component needs the complete state to initialize, for instance to push the initial sediment onto the sediment buffer, it implements `initialize!(::Input, input, state)`, which is called after all parts are created. A model's `initial_state` is then a one-liner:

```julia
initial_state(input::AbstractInput) = Common.initial_state(State, input)
```

## Contents

```@contents
Pages = ["tag.md", "boxes.md", "time.md", "facies.md",
         "production.md", "cellular-automata.md", "waterdepth.md", "hdf5.md",
         "sediment_buffer.md"]
Depth = 1
```

## Common Definitions

``` {.julia file=src/Components/Common.jl}
module Common
export @u_str, Quantity, Amount, Time, Location, Rate, Intensity, Height, Sediment
export AbstractFacies, AbstractInput, AbstractState, AbstractOutput, set_attribute
export InputPart, StatePart, FaciesPart, write_header
export OutputSpec
export Box, box_axes, Boundary, Coast, Periodic, Reflected, TimeProperties
export in_units_of
export Model
export record_type
export Size
export n_steps
export AbstractProduction, NoProduction, BenthicProduction, PelagicProduction, production_profile, InterpolatedProduction

using ModuleMixins.Procedural: record_type, parts, map_components
using Unitful
using StaticArrays

using ...CarboKitten: Model, AbstractInput, AbstractState, AbstractFacies, AbstractOutput, OutputSpec, set_attribute, n_steps
using ...CarboKitten: InputPart, StatePart, FaciesPart
using ...BoundaryTrait
using ...Config: TimeProperties
using ...Boxes: Box, box_axes
using ...Utility: in_units_of
using ...Production

const Amount = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")
const Height = typeof(1.0u"m")
const Location = typeof(1.0u"m")
const Rate = typeof(1.0u"m/Myr")
const Intensity = typeof(1.0u"W/m^2")
const Sediment = typeof(1.0u"m")

"""
    write_header(part, input, output)

Write the attributes of a single component to the output header. Components
extend this method for their own `Input` part.
"""
write_header(part, input::AbstractInput, output::AbstractOutput) = nothing

"""
    write_header(input, output)

Write the header for all components of a model.
"""
write_header(input::AbstractInput, output::AbstractOutput) =
    foreach(p -> write_header(p, input, output), parts(input))

"""
    initial_state(part, input)

The initial `State` part of a single component, or `nothing` if the component
has no state. Components extend this method for their own `Input` part.
"""
initial_state(part::InputPart, input::AbstractInput) = nothing

"""
    initialize!(part, input, state)

Finish initializing the `state`, once all its parts are created. Components
that need the complete state to initialize, extend this method for their own
`Input` part.
"""
initialize!(part::InputPart, input::AbstractInput, state::AbstractState) = nothing

"""
    initial_state(S, input)

The initial state of type `S` for a model. The parts of the state are
collected from all components, after which `initialize!` is called for each
component, in dependency order.
"""
function initial_state(::Type{S}, input::AbstractInput) where {S <: AbstractState}
    # `map_components` creates a `Composite` without knowledge of roles.
    # Converting to `S` makes role lookups type stable.
    state = S(parts(map_components(p -> initial_state(p, input), input)))
    foreach(p -> initialize!(p, input, state), parts(input))
    return state
end

end
```

``` {.julia file=src/Components.jl}
module Components

export Tag, TimeIntegration, Boxes, WaterDepth, FaciesBase, Production,
    CAProduction, CAFeedback, CellularAutomaton, H5Writer, ActiveLayer, SedimentBuffer,
    CopyingSedimentBuffer, CircularSedimentBuffer,
    ActiveLayerOnshore, Denudation, InitialSediment, Output, Diagnostics

using ModuleMixins.Syntactic: @component

include("Components/Common.jl")
include("Components/Tag.jl")
include("Components/Diagnostics.jl")
include("Components/TimeIntegration.jl")
include("Components/Boxes.jl")
include("Components/WaterDepth.jl")
include("Components/FaciesBase.jl")
include("Components/Production.jl")
include("Components/CellularAutomaton.jl")
include("Components/CAProduction.jl")
include("Components/CAFeedback.jl")

include("Components/SedimentBuffer.jl")
include("Components/CopyingSedimentBuffer.jl")
include("Components/CircularSedimentBuffer.jl")
include("Components/InitialSediment.jl")
include("Components/ActiveLayer.jl")
include("Components/Denudation.jl")

include("Components/Output.jl")

list_components() = filter(
    c -> c isa Module && c !== Components && c !== Common,
    map(eval, names(Components)))

end
```
