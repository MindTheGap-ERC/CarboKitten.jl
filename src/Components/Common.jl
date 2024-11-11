# ~/~ begin <<docs/src/components/components.md#src/Components/Common.jl>>[init]
module Common
export @u_str, Quantity, Amount, Time, Location, Rate, Intensity, Height
export AbstractFacies, AbstractInput, AbstractState, AbstractFrame
export Box, axes, Boundary, Shelf, Periodic, Reflected, TimeProperties
export in_units_of
export Model
export @for_each

using ModuleMixins
using Unitful
using CarboKitten.BoundaryTrait
using CarboKitten.Config: TimeProperties
using CarboKitten.Boxes: Box, axes
using CarboKitten.Utility: in_units_of

const Amount = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")
const Height = typeof(1.0u"m")
const Location = typeof(1.0u"m")
const Rate = typeof(1.0u"m/Myr")
const Intensity = typeof(1.0u"W/m^2")

abstract type AbstractFacies end
abstract type AbstractInput end
abstract type AbstractState end
abstract type AbstractFrame end

struct Model{M} end

end
# ~/~ end
