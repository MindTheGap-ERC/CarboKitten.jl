module CarboKitten

include("./BoundaryTrait.jl")
include("./Vectors.jl")
include("./Config.jl")
include("./Stencil.jl")
include("./SedimentStack.jl")
include("./Burgess2013.jl")
include("./Utility.jl")
# include("./BS92.jl")
include("./CaProd.jl")
include("./Visualization.jl")

module Model
include("./Model/ALCAPS.jl")
end

end # module CarboKitten
