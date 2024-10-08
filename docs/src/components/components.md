# Model Components

Each model in CarboKitten is composed out of elementary parts. These form a hierarchy of modules that inherit from each other.

```mermaid
graph TD;
    classDef nyi fill:#888;

    time[Time Integration]
    facies[Facies Base]
    box[Box]
    water_depth[Water Depth]
    uniform_production[Production]
    cellular_automaton[Cellular Automaton]
    ca_production[CA Production]

    sediment_buffer[Sediment Buffer]
    active_layer_transport[Active Layer Transport]
    denudation[Denudation]
    class sediment_buffer,active_layer_transport,denudation nyi;

    time --> water_depth
    box --> water_depth

    water_depth --> uniform_production
    facies --> uniform_production

    box --> cellular_automaton
    facies --> cellular_automaton

    cellular_automaton --> ca_production
    uniform_production --> ca_production

    box --> sediment_buffer
    sediment_buffer --> active_layer_transport
    water_depth --> active_layer_transport
    sediment_buffer --> denudation
    water_depth --> denudation
```

## Contents

```@contents
Pages = ["boxes.md", "time.md", "facies.md",
         "production.md", "cellular-automata.md", "waterdepth.md"]
Depth = 1
```

## Common Definitions

``` {.julia file=src/Components/Common.jl}
module Common
    export @u_str, Amount, Time, Location, Rate, Intensity, Height
    export AbstractFacies, AbstractInput, AbstractState
    export Box, axes, Boundary, Shelf, Periodic, Reflected, TimeProperties

    using Unitful
    using CarboKitten.BoundaryTrait
    using CarboKitten.Config

    const Amount = typeof(1.0u"m")
    const Time = typeof(1.0u"Myr")
    const Height = typeof(1.0u"m")
    const Location = typeof(1.0u"m")
    const Rate = typeof(1.0u"m/Myr")
    const Intensity = typeof(1.0u"W/m^2")

    abstract type AbstractFacies end
    abstract type AbstractInput end
    abstract type AbstractState end
end
```

``` {.julia file=src/Components.jl}
module Components

export TimeIntegration, Boxes, WaterDepth, FaciesBase, UniformProduction, CAProduction, CellularAutomaton

using ModuleMixins: @compose

include("Components/Common.jl")
include("Components/TimeIntegration.jl")
include("Components/Boxes.jl")
include("Components/WaterDepth.jl")
include("Components/FaciesBase.jl")
include("Components/Production.jl")
include("Components/CellularAutomaton.jl")
include("Components/CAProduction.jl")

end
```
