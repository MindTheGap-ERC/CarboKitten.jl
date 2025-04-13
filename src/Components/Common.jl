# ~/~ begin <<docs/src/components/components.md#src/Components/Common.jl>>[init]
module Common
export @u_str, Quantity, Amount, Time, Location, Rate, Intensity, Height
export AbstractFacies, AbstractInput, AbstractState, AbstractFrame
export Box, box_axes, Boundary, Coast, Periodic, Reflected, TimeProperties
export in_units_of
export Model
export @for_each
export Size
export Frame

using ModuleMixins
using Unitful
using StaticArrays

using ...CarboKitten: Model
using ...BoundaryTrait
using ...Config: TimeProperties
using ...Boxes: Box, box_axes
using ...Utility: in_units_of

const Amount = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")
const Height = typeof(1.0u"m")
const Location = typeof(1.0u"m")
const Rate = typeof(1.0u"m/Myr")
const Intensity = typeof(1.0u"W/m^2")
const Sediment = typeof(1.0u"m")

abstract type AbstractFacies end
abstract type AbstractInput end
abstract type AbstractState end
abstract type AbstractFrame end

@kwdef struct Frame
    disintegration::Union{Array{Sediment,3},Nothing} = nothing   # facies, x, y
    production::Union{Array{Sediment,3},Nothing} = nothing
    deposition::Union{Array{Sediment,3},Nothing} = nothing
end

end
# ~/~ end
