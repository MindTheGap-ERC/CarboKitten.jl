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
include("./Visualization.jl")

module Transport
include("./Transport/ActiveLayer.jl")
end

module Model
include("./Model/ALCAPS.jl")
include("./Model/WithDenudation.jl")
end

end # module CarboKitten
