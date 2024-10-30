# Utility functions

## Select Iterator
In many cases our model is updating some state as an iterator. I want to be able to select arbitrary slices from that iterator. The `Select` object contains the main iterable and a selection iterable that should yield integers. When iterated upon, we repeatedly use `Iterators.dropwhile` to get to the next index.

``` {.julia file=src/Utility.jl}
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

<<utility>>

end
```

## Range finder

The `RangeFinder` iterator is used in our algorithm to trace hiatus in the sediment history. This iterator consumes an iterator of booleans and yields values of type `UnitRange`, giving all ranges for which the input sequence is true consecutively.

``` {.julia #utility-spec}
a = [false, true, true, false, true]
@test collect(find_ranges(a)) == [2:3, 5:5]
```

``` {.julia #utility}
struct RangeFinder
	v::AbstractVector{Bool}
end

Base.iterate(r::RangeFinder) = iterate(r, 1)

function Base.iterate(r::RangeFinder, i::Union{Int, Nothing})
	isnothing(i) && return nothing
	a = findnext(r.v, i)
	isnothing(a) && return nothing
	b = findnext(!, r.v, a)
	isnothing(b) && return (a:length(r.v)), nothing
	return (a:b-1), b
end

Base.eltype(r::RangeFinder) = UnitRange{Int}
Base.IteratorSize(::Type{RangeFinder}) = Base.SizeUnknown()

"""
    find_ranges(v::AbstractVector{Bool})

Take a vector of bools, returns an iterator over all ranges for
which the vector is `true`.
"""
find_ranges(v::AbstractVector{Bool}) = RangeFinder(v)

export find_ranges
```

## Tagging a sequence of vectors

The `enumerate_seq` iterator is also used in the algorithm to trace hiatus in sediment accumulation. Here the task is to enumerate a nested sequence while preserving the nested structure.

``` {.julia #utility-spec}
a = [[:a, :b], [:c]]
@test collect(enumerate_seq(a)) == [[(1, :a), (2, :b)], [(3, :c)]]
```

``` {.julia #utility}
struct TagVectors{T}
	vectors::T
end

function Base.iterate(tv::TagVectors{T}) where T
	tag = 1
	x = iterate(tv.vectors)
	isnothing(x) && return nothing
	(v_, nit) = x
	v = collect(v_)
	n = length(v)
	return zip(tag:tag+n-1, v) |> collect, (tag+n, nit)
end

function Base.iterate(tv::TagVectors{T}, st::Tuple{Int,U}) where {T, U}
	(tag, it) = st
	isnothing(it) && return nothing
	x = iterate(tv.vectors, it)
	isnothing(x) && return nothing
	(v_, nit) = x
	v = collect(v_)
	n = length(v)
	return zip(tag:tag+n-1, v) |> collect, (tag+n, nit)
end

Base.IteratorSize(::Type{TagVectors{T}}) where T = Base.IteratorSize(T)
Base.size(r::TagVectors{T}) where T = size(r.vectors)
Base.length(r::TagVectors{T}) where T = length(r.vectors)
Base.eltype(::Type{TagVectors{T}}) where T = Any
# This should be:
# Iterators.Zip{Tuple{UnitRange{Int},Vector{eltype(eltype(T))}}}
# But that doesn't work

"""
    enumerate_seq(s)

Enumerates an iterator of sequences, such that the following equivalence
holds:

    enumerate(flatten(s)) == flatten(enumerate_seq(s))
"""
enumerate_seq(s::T) where T = TagVectors{T}(s)

export enumerate_seq
```

``` {.julia file=test/UtilitySpec.jl}
@testset "CarboKitten.Utility" begin
    using CarboKitten.Utility

    <<utility-spec>>
end
```

## Reading data files

``` {.julia file=src/DataSets.jl}
module DataSets

export read_tsv, read_csv

using DelimitedFiles: readdlm
using DataFrames
using CategoricalArrays
using CSV
using Pkg.Artifacts
using Unitful

function read_tsv(filename)
    data, header = readdlm(filename, '\t', header=true)
    return DataFrame(data, vec(header))
end

function read_csv(filename)
    return DataFrame(CSV.File(filename))
end

function artifact_dir()
    subfolder = first(readdir(artifact"data"))
    return joinpath(artifact"data", subfolder)
end

function bosscher_schlager_1992()
    dir = artifact_dir()
    filename = joinpath(dir, "Bosscher1992", "bs92-sealevel-curve.csv")
    df = read_csv(filename)
    return DataFrame(time=df.time * u"yr", sealevel=df.depth * u"m")
end

function miller_2020()
    dir = artifact_dir()
    filename = joinpath(dir, "Miller2020", "Cenozoic_sea_level_reconstruction.tab")
    df = read_tsv(filename)
    return DataFrame(
        time=-df[!,4] * u"kyr",
        sealevel=df[!,7] * u"m",
        reference=categorical(df[!,2]))
end

end
```
