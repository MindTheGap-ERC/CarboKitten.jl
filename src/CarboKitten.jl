module CarboKitten

include("./DSL.jl")
include("./BoundaryTrait.jl")
include("./Vectors.jl")
include("./Config.jl")
include("./Stencil.jl")
include("./SedimentStack.jl")
include("./Utility.jl")
include("./Burgess2013.jl")

include("./Denudation.jl")
include("./CaProd.jl")
include("./Visualization.jl")

module Transport
include("./Transport/ActiveLayer.jl")
end

module Model
using ModuleMixins: @compose
using CarboKitten.Components.Common
using CarboKitten.Components

include("./Model/BS92.jl")
include("./Model/ALCAPS.jl")
include("./Model/WithDenudation.jl")
end

include("./Components.jl")

end # module CarboKitten
