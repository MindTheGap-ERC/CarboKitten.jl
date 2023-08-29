# ~/~ begin <<docs/src/utility.md#src/Utility.jl>>[init]
module Utility

export select

struct Select
    iter
    selection
end

function select(it, sel)
    Select(enumerate(it), sel)
end

function Base.iterate(s::Select)
    x = iterate(s.selection)
    if x !== nothing
        (idx, selstate) = x
        ((_, value), rest) = Iterators.peel(Iterators.dropwhile(((i, y),) -> i != idx, s.iter))
        return (value, (selstate, rest))
    else
        return nothing
    end
end

function Base.iterate(s::Select, state)
    (selstate, rest) = state
    x = iterate(s.selection, selstate)
    if x !== nothing
        (idx, selstate) = x
        ((_, value), rest) = Iterators.peel(Iterators.dropwhile(((i, y),) -> i != idx, s.iter))
        return (value, (selstate, rest))
    else
        return nothing
    end
end

end
# ~/~ end