module CarboKitten

include("./BoundaryTrait.jl")
include("./Vectors.jl")
include("./Config.jl")
include("./Stencil.jl")
include("./SedimentStack.jl")
include("./Burgess2013.jl")
include("./Utility.jl")
#include("./BS92.jl")
include("./InputConfig.jl")
include("./Denudation.jl")

include("./CaProd.jl")
include("./CaProdErosion.jl")

#include("./Visualization.jl")


include("./Visualization.jl")

module Transport
include("./Transport/ActiveLayer.jl")
end

module Model
include("./Model/ALCAPS.jl")
end

end # module CarboKitten
