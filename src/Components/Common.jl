# ~/~ begin <<docs/src/components/components.md#src/Components/Common.jl>>[init]
module Common
export @u_str, Quantity, Amount, Time, Location, Rate, Intensity, Height, Sediment
export AbstractFacies, AbstractInput, AbstractState, AbstractOutput, set_attribute
export OutputSpec
export Box, box_axes, Boundary, Coast, Periodic, Reflected, TimeProperties
export in_units_of
export Model
export @for_each
export Size
export n_steps

using ModuleMixins
using Unitful
using StaticArrays

using ...CarboKitten: Model, AbstractInput, AbstractState, AbstractOutput, OutputSpec, set_attribute, n_steps
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

end
# ~/~ end
