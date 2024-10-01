# Model Components

Each model in CarboKitten is composed out of elementary parts.

``` {.julia file=src/Components/Common.jl}
module Common
    export @u_str, Amount, Time, Location, Rate, Intensity, Height
    export AbstractFacies, AbstractInput, AbstractState
    export Box, axes, Boundary, Shelf, Periodic, Reflected

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

1. [Time](./time.md)
2. [Water depth](./waterdepth.md)
4. [Facies](./facies.md)
3. [Production](./production.md)
4. [Cellular Automata](./cellular-automata.md)

``` {.julia file=src/Components.jl}
include("Components/Common.jl")
include("Components/TimeIntegration.jl")
include("Components/Boxes.jl")
include("Components/WaterDepth.jl")
include("Components/FaciesBase.jl")
include("Components/Production.jl")
include("Components/CellularAutomata.jl")
```

