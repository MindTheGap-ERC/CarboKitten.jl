# Range finder

The `RangeFinder` iterator is used in our algorithm to trace hiatus in the sediment history. This iterator consumes an iterator of booleans and yields values of type `UnitRange`, giving all ranges for which the input sequence is true consecutively.

``` {.julia #algoritms-spec}
a = [false, true, true, false, true]
@test collect(find_ranges(a)) == [2:3, 5:5]
```

``` {.julia file=src/Algorithms/RangeFinder.jl}
module RangeFinder

struct _RangeFinder{T}
	v::T
end

Base.iterate(r::_RangeFinder{T}) where {T} = iterate(r, 1)

function Base.iterate(r::_RangeFinder{T}, i::Union{Int, Nothing}) where {T}
	isnothing(i) && return nothing
	a = findnext(r.v, i)
	isnothing(a) && return nothing
	b = findnext(!, r.v, a)
	isnothing(b) && return (a:length(r.v)), nothing
	return (a:b-1), b
end

Base.eltype(r::_RangeFinder{T}) where {T} = UnitRange{Int}
Base.IteratorSize(::Type{_RangeFinder{T}}) where {T} = Base.SizeUnknown()

"""
    find_ranges(v::AbstractVector{Bool})

Take a vector of bools, returns an iterator over all ranges for
which the vector is `true`.
"""
find_ranges(v::T) where {T <: AbstractVector{Bool}} = _RangeFinder{T}(v)

export find_ranges

end
```
