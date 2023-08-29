---
title: CarboKitten
subtitle: utility functions
---

# Utility functions

## Select Iterator
In many cases our model is updating some state as an iterator. I want to be able to select arbitrary slices from that iterator. The `Select` object contains the main iterable and a selection iterable that should yield integers. When iterated upon, we repeatedly use `Iterators.dropwhile` to get to the next index.

``` {.julia file=src/Utility.jl}
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
```
