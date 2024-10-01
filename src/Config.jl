# ~/~ begin <<docs/src/boxes.md#src/Config.jl>>[init]
module Config

export AbstractBox, Box, TimeProperties, axes

using ..BoundaryTrait
using ..Vectors

using Unitful
using Unitful.DefaultSymbols

# ~/~ begin <<docs/src/boxes.md#config-types>>[init]
abstract type AbstractBox{BT} end

struct Box{BT} <: AbstractBox{BT}
    grid_size::NTuple{2,Int}
    phys_scale::typeof(1.0m)
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
# ~/~ begin <<docs/src/boxes.md#config-types>>[1]
abstract type AbstractTimeProperties end

@kwdef struct TimeProperties <: AbstractTimeProperties
    Δt::typeof(1.0u"Myr")
    steps::Int
    write_interval::Int
end
# ~/~ end

end
# ~/~ end
