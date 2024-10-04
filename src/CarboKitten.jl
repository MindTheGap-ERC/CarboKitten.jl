module CarboKitten

include("./BoundaryTrait.jl")
include("./Vectors.jl")
include("./Config.jl")
include("./Stencil.jl")
include("./SedimentStack.jl")
include("./Utility.jl")
include("./Burgess2013.jl")

include("./Denudation.jl")
include("./CaProd.jl")

module Transport
include("./Transport/ActiveLayer.jl")
end

include("./Components.jl")

module Model
using ModuleMixins: @compose
using CarboKitten.Components.Common
using CarboKitten.Components: Production, WaterDepth, TimeIntegration, Boxes, FaciesBase

include("./Model/BS92.jl")
include("./Model/ALCAPS.jl")
include("./Model/WithDenudation.jl")
end

include("./Export.jl")
include("./Visualization.jl")

end # module CarboKitten
