module CarboKitten

using TerminalLoggers: TerminalLogger
using Logging

function init()
    global_logger(TerminalLogger(right_justify=80))
    @info """
    # Welcome to CarboKitten!
    """
end

include("./BoundaryTrait.jl")
include("./Vectors.jl")
include("./Config.jl")
include("./Boxes.jl")
include("./Stencil.jl")
include("./SedimentStack.jl")
include("./Utility.jl")
include("./Skeleton.jl")

include("./Burgess2013.jl")

include("./Denudation.jl")

module Transport
include("./Transport/ActiveLayer.jl")
end

include("./Components.jl")

module Model
using ModuleMixins: @compose
using CarboKitten.Components.Common
using CarboKitten.Components

include("./Model/BS92.jl")
include("./Model/CAP.jl")
include("./Model/ALCAP2.jl")
include("./Model/WithDenudation.jl")
end

include("./Export.jl")
include("./Visualization.jl")

const run = Components.H5Writer.run

end # module CarboKitten
