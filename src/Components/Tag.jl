# ~/~ begin <<docs/src/components/tag.md#src/Components/Tag.jl>>[init]
@compose module Tag
using ..Common
using HDF5

@kwdef struct Input <: AbstractInput
    tag::String = "untagged run"
end

function write_header(fid, input::AbstractInput)
    attr = attributes(fid["input"])
    attr["tag"] = input.tag
end
end
# ~/~ end