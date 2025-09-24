# ~/~ begin <<docs/src/components/boxes.md#src/Components/Boxes.jl>>[init]
@compose module Boxes
using ..Common

@kwdef struct Input <: AbstractInput
    box::Box
end

function write_header(input::AbstractInput, output::AbstractOutput)
    x, y = box_axes(input.box)
    set_attribute(output, "x", collect(x) |> in_units_of(u"m"))
    set_attribute(output, "y", collect(y) |> in_units_of(u"m"))
    set_attribute(output, "grid_size_x", input.box.grid_size[1])
    set_attribute(output, "grid_size_y", input.box.grid_size[2])
    set_attribute(output, "phys_scale", input.box.phys_scale |> in_units_of(u"m"))
end
end
# ~/~ end
