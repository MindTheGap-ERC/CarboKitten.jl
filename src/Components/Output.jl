# ~/~ begin <<docs/src/components/output.md#src/Components/Output.jl>>[init]
@component module Output
using ..Common

# ~/~ begin <<docs/src/components/output.md#output-spec>>[init]
@kwdef struct Input <: InputPart
    output = Dict(:full => OutputSpec((:, :), 1))
end
# ~/~ end
end
# ~/~ end
