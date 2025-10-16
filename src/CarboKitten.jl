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
function get_logger end

abstract type AbstractBox{BT} end
abstract type AbstractInput end
abstract type AbstractOutput end
abstract type AbstractState end

const Slice2 = NTuple{2,Union{Int,Colon,UnitRange{Int}}}

@kwdef struct OutputSpec
    slice::Slice2 = (:, :)
    write_interval::Int = 1
end

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
    y_axis = (0:(box.grid_size[2]-1)) .* box.phys_scale
    x_axis = (0:(box.grid_size[1]-1)) .* box.phys_scale
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
    n_steps(input)

Returns the number of steps from a given input.
"""
function n_steps end


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

include("./RunModel.jl")

include("./Denudation.jl")

module Transport
include("./Transport/ActiveLayer.jl")
include("./Transport/DifferentialOperators.jl")
include("./Transport/Solvers.jl")
include("./Transport/Advection.jl")
end

function set_attribute end

include("./Components.jl")

module Output

include("./Output/Abstract.jl")
include("./Output/RunModel.jl")
include("./Output/H5Writer.jl")
include("./Output/MemoryWriter.jl")

using .Abstract: Frame, frame_writer, state_writer
export Frame, frame_writer, state_writer

end

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
include("./Visualization.jl")
include("./Testing.jl")

using .Components.Common: in_units_of, @u_str
using .Output.Abstract: OutputSpec, new_output
using .Output.MemoryWriter: MemoryOutput
using .Models: BS92, CAP, ALCAP
using .BoundaryTrait: Boundary, Coast, Periodic, Reflected
using GeometryBasics: Vec2

export run_model, Box, box_axes, TimeProperties, time_axis,
    Model, BS92, CAP, ALCAP, in_units_of, @u_str,
    AbstractBox, Boundary, Coast, Periodic, Reflected,
    Vec2, OutputSpec, MemoryOutput, new_output, n_steps

end # module CarboKitten
