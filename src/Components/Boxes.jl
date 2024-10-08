# ~/~ begin <<docs/src/components/boxes.md#src/Components/Boxes.jl>>[init]
@compose module Boxes
using ..Common

@kwdef struct Input <: AbstractInput
    box::Box
end

function write_header(fid, input::AbstractInput)
    x, y = Common.axes(input.box)

    gid = fid["input"]
    gid["x"] = collect(x) |> in_units_of(u"m")
    gid["y"] = collect(y) |> in_units_of(u"m")
end
end
# ~/~ end
