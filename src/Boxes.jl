# ~/~ begin <<docs/src/components/boxes.md#src/Boxes.jl>>[init]
module Boxes

using ..BoundaryTrait
using ..Vectors
using Unitful
using Unitful.DefaultSymbols

export AbstractBox, Box, axes

# ~/~ begin <<docs/src/components/boxes.md#box-type>>[init]
abstract type AbstractBox{BT} end

struct Box{BT} <: AbstractBox{BT}
    grid_size::NTuple{2,Int}
    phys_scale::typeof(1.0u"m")
    phys_size::Vec2

    function Box{BT}(;grid_size::NTuple{2, Int}, phys_scale::Quantity{Float64, 𝐋, U}) where {BT <: Boundary{2}, U}
        new{BT}(grid_size, phys_scale, phys_size(grid_size, phys_scale))
    end
end

function axes(box::Box)
	y_axis = (0:(box.grid_size[2] - 1)) .* box.phys_scale
	x_axis = (0:(box.grid_size[1] - 1)) .* box.phys_scale
	return x_axis, y_axis
end

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