# ~/~ begin <<docs/src/components/boxes.md#src/Boxes.jl>>[init]
module Boxes

import ..CarboKitten
import ..BoundaryTrait

Base.@deprecate_binding AbstractBox CarboKitten.AbstractBox
Base.@deprecate_binding Box CarboKitten.Box
Base.@deprecate_binding box_axes CarboKitten.box_axes
Base.@deprecate_binding axes CarboKitten.box_axes false

Base.@deprecate_binding Boundary BoundaryTrait.Boundary false
Base.@deprecate_binding Reflected BoundaryTrait.Reflected false
Base.@deprecate_binding Periodic BoundaryTrait.Periodic false
Base.@deprecate_binding Constant BoundaryTrait.Constant false
Base.@deprecate_binding Coast BoundaryTrait.Coast false
Base.@deprecate_binding Shelf BoundaryTrait.Coast false

end
# ~/~ end
