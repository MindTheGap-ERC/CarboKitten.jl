# ~/~ begin <<docs/src/components/tag.md#src/Components/Tag.jl>>[init]
@compose module Tag
using ..Common
using HDF5

@kwdef struct Input <: AbstractInput
    tag::String = "untagged run"
end

function write_header(input::AbstractInput, output::AbstractOutput)
    set_attribute(output, "tag", input.tag)
end
end
# ~/~ end
