module CarboKitten

using TerminalLoggers: TerminalLogger
using Logging
using Unitful

function init()
    global_logger(TerminalLogger(right_justify=80))
    @info """
    # Welcome to CarboKitten!
    """
end

struct Model{M} end
function run_model end

abstract type AbstractBox{BT} end

"""
    Box{BoundaryType}(grid_size, phys_scale)

Stores the spatial properties of the simulation grid, where `grid_size` is a tuple of two integers and `phys_scale` is a `Quantity` of dimension `Length`.
"""
@kwdef struct Box{BT} <: AbstractBox{BT}
    grid_size::NTuple{2,Int}
    phys_scale::typeof(1.0u"m")
end

"""
    box_axes(box::Box)

Return the `x` and `y` components of the box grid in meters. The grid starts at
`0.0u"m"` and runs to `(box.grid_size .- 1) .* box.phys_scale`.
"""
function box_axes(box::Box)
	y_axis = (0:(box.grid_size[2] - 1)) .* box.phys_scale
	x_axis = (0:(box.grid_size[1] - 1)) .* box.phys_scale
	return x_axis, y_axis
end

abstract type AbstractTimeProperties end

"""
    TimeProperties(t0, Δt, steps, write_interval)

Stores properties of the time integration. Here, `t0` and `Δt` should be a
`Quantity`, `steps` is the number of integration steps, and `write_interval` the
number of steps between writing output.
"""
@kwdef struct TimeProperties <: AbstractTimeProperties
    t0::typeof(1.0u"Myr") = 0.0u"Myr"
    Δt::typeof(1.0u"Myr")
    steps::Int
end

n_writes(time::TimeProperties) = time.steps

"""
    time_axis(time::TimeProperties)

Retrieve the time values for which output was/will be written. Returns a range.
"""
time_axis(time::TimeProperties) = (0:n_writes(time)) .* time.Δt .+ time.t0

include("./BoundaryTrait.jl")
include("./Vectors.jl")
include("./Boxes.jl")
include("./Config.jl")
include("./Stencil.jl")
include("./SedimentStack.jl")
include("./Utility.jl")
include("./DataSets.jl")
include("./Skeleton.jl")
include("./OutputData.jl")

include("./Denudation.jl")

module Transport
include("./Transport/ActiveLayer.jl")
include("./Transport/DifferentialOperators.jl")
include("./Transport/Solvers.jl")
include("./Transport/Advection.jl")
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
include("./Models/WithoutCA.jl")

end

include("./Export.jl")
include("./MemoryWriter.jl")
include("./Visualization.jl")
include("./Testing.jl")

using .Components.Common: in_units_of, @u_str
using .OutputData: OutputSpec, new_output
using .MemoryWriter: MemoryOutput
using .Models: BS92, CAP, ALCAP
using .BoundaryTrait: Boundary, Coast, Periodic, Reflected
using GeometryBasics: Vec2

export run_model, Box, box_axes, TimeProperties, time_axis,
       Model, BS92, CAP, ALCAP, in_units_of, @u_str,
       AbstractBox, Boundary, Coast, Periodic, Reflected,
       Vec2, OutputSpec, MemoryOutput, new_output

end # module CarboKitten
