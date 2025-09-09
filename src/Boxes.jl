# ~/~ begin <<docs/src/components/boxes.md#src/Boxes.jl>>[init]
module Boxes

using ..BoundaryTrait
using ..CarboKitten: AbstractBox, Box, box_axes
using ..Vectors

using Unitful
using Unitful.DefaultSymbols

export AbstractBox, Box, box_axes

# ~/~ begin <<docs/src/components/boxes.md#box-type>>[init]
const axes = box_axes

phys_size(grid_size, phys_scale) = (
    x = grid_size[1] * (phys_scale / m |> NoUnits),
    y = grid_size[2] * (phys_scale / m |> NoUnits))
# ~/~ end
# ~/~ begin <<docs/src/components/boxes.md#vector-offset>>[init]
Base.in(a::Vec2, box::Box) =
    a.x >= 0.0 && a.x < box.phys_size.x && a.y >= 0.0 && a.y < box.phys_size.y

function offset(box::AbstractBox{Reflected{2}}, a::Vec2, Δa::Vec2)
    clip(i, a, b) = (i < a ? a + a - i : (i > b ? b + b - i : i))
    (x=clip(a.x+Δa.x, 0.0, box.phys_size.x)
    ,y=clip(a.y+Δa.y, 0.0, box.phys_size.y))
end

function offset(box::AbstractBox{Periodic{2}}, a::Vec2, Δa::Vec2)
    (x=mod(a.x+Δa.x, box.phys_size.x)
    ,y=mod(a.y+Δa.y, box.phys_size.y))
end

function offset(box::AbstractBox{Constant{2,Value}}, a::Vec2, Δa::Vec2) where Value
    b = a + Δa
    if b ∉ box
        nothing
    else
        b
    end
end

function offset(box::AbstractBox{Shelf}, a::Vec2, Δa::Vec2)
    b = a + Δa
    if b.x < 0.0 || b.x >= box.phys_size.x
        nothing
    else
        (x=b.x, y=mod(b.y, box.phys_size.y))
    end
end
# ~/~ end

end
# ~/~ end