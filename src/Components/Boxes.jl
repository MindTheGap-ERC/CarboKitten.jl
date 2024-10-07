# ~/~ begin <<docs/src/components/boxes.md#src/Components/Boxes.jl>>[init]
@compose module Boxes
    using ..Common

    @kwdef struct Input <: AbstractInput
        box::Box
    end
end
# ~/~ end
