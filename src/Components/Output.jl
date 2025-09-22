# ~/~ begin <<docs/src/components/output.md#src/Components/Output.jl>>[init]
@compose module Output
    # ~/~ begin <<docs/src/components/output.md#output-spec>>[init]
    @kwdef struct Input <: AbstractInput
        output = Dict(:full => OutputSpec((:,:), 1))
    end
    # ~/~ end
end
# ~/~ end