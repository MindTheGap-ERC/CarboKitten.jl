# ~/~ begin <<docs/src/components/components.md#src/Components.jl>>[init]
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
# ~/~ end
