# ~/~ begin <<docs/src/algorithms/stratigraphic_column.md#src/Algorithms/StratigraphicColumn.jl>>[init]
module StratigraphicColumn

export stratigraphic_column!

"""
    stratigraphic_column!(a)

Computes the stratigraphic column, given a time series of sediment
events. The input sequence will contain positive numbers for deposits
and negative numbers for disintegrates. The column is computed in-place.

Returns the mutated input.
"""
function stratigraphic_column!(a::AbstractVector{T}) where {T}
    for i in 2:length(a)
        a[i] += a[i-1]
    end
    for i in length(a):-1:2
        a[i-1] = min(a[i-1], a[i])
        a[i] -= a[i-1]
    end
    return a
end

end
# ~/~ end
