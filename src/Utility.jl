# ~/~ begin <<docs/src/utility.md#src/Utility.jl>>[init]
module Utility

export select, in_units_of
using Unitful

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

function in_units_of(unit)
    function magnitude(a::AbstractArray{Quantity{RT, NoDims, U}, dim}) where {RT <: Real, U, dim}
        return a .|> NoUnits
    end

    function magnitude(a::AbstractArray{RT, dim}) where {RT <: Real, dim}
        return a
    end

    function magnitude(a::RT) where {RT <: Real}
        return a
    end

    function magnitude(a::Quantity{RT, NoDims, U}) where {RT <: Real, U}
        return a |> NoUnits
    end

    function (x)
        x / unit |> magnitude
    end
end

end
# ~/~ end