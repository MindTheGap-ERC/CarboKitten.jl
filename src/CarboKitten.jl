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
include("./Boxes.jl")
include("./Config.jl")
include("./Stencil.jl")
include("./SedimentStack.jl")
include("./Utility.jl")
include("./DataSets.jl")
include("./Skeleton.jl")

include("./Burgess2013.jl")

include("./Denudation.jl")

module Transport
include("./Transport/ActiveLayer.jl")
end

include("./Components.jl")

module Models
using ModuleMixins: @compose
using CarboKitten.Components.Common
using CarboKitten.Components

include("./Models/BS92.jl")
include("./Models/CAP.jl")
include("./Models/ALCAP.jl")
include("./Models/WithDenudation.jl")
end

include("./Export.jl")
include("./Visualization.jl")

using .Components.H5Writer: run_model
using .Boxes: Box, box_axes
using .Config: TimeProperties
using .Components.TimeIntegration: TimeProperties, time_axis
using .Components.Common: Model, in_units_of, @u_str
using .Models: BS92, CAP, ALCAP
using .BoundaryTrait: Boundary, Coast, Periodic, Reflected

export run_model, Box, box_axes, TimeProperties, time_axis,
       Model, BS92, CAP, ALCAP, in_units_of, @u_str,
       Boundary, Coast, Periodic, Reflected

end # module CarboKitten
