# ~/~ begin <<docs/src/boxes.md#src/Config.jl>>[init]
module Config

export TimeProperties

using ..BoundaryTrait
using ..Vectors

using ..Boxes: Box, axes
export Box, axes

using Unitful
using Unitful.DefaultSymbols

# ~/~ begin <<docs/src/boxes.md#config-types>>[init]
abstract type AbstractTimeProperties end

@kwdef struct TimeProperties <: AbstractTimeProperties
    t0::typeof(1.0u"Myr") = 0.0u"Myr"
    Δt::typeof(1.0u"Myr")
    steps::Int
    write_interval::Int = 1
end
# ~/~ end

end
# ~/~ end
