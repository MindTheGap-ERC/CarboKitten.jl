# ~/~ begin <<docs/src/stencils.md#src/Stencil.jl>>[init]
module Stencil

using StaticArrays

using ..Boxes: AbstractBox
using ..BoundaryTrait

export stencil, convolution

# ~/~ begin <<docs/src/stencils.md#stencil-operation>>[init]
"""
    stencil!(f, Boundary, Size, out, inp...)

Performs stencil operation. The function `f` should take a number of
abstract arrays (same as the number of `inp` values),
and return a single value of the same type as elements in `out`.
The `Size` parameter should be a type parameter size of the stencil,
e.g. `Size(3, 3)` to get a 3 by 3 stencil.

Prefer to use this version over the older implementations.
"""
function stencil!(f::F, ::Type{BT}, ::Size{sz}, out, inp...) where {F, dim, sz, BT <: Boundary{dim}}
    stencil_multi!(f, BT, Size(sz), out, inp)
end

function stencil_multi!(f::F, ::Type{BT}, ::Size{sz}, out, inp) where {F, dim, sz, BT <: Boundary{dim}}
    @assert(
        all(size(a) == size(out) for a in inp),
        "inputs have wrong shape: $([size(a) for a in inp]), should be $(size(out))")

    center = CartesianIndex((div.(sz, 2) .+ 1)...)
    for i in eachindex(IndexCartesian(), out)
        nb = (SArray{Tuple{sz...}}(
                get_bounded(BT, a, i + j - center)
                for j in CartesianIndices(sz))
              for a in inp)
        out[i] = f(nb...)
    end
    return out
end


function stencil(::Type{TIn}, ::Type{TOut}, ::Type{BT}, n::NTuple{dim,Int}, f::Function) where {TIn, TOut, dim, BT <: Boundary{dim}}
    m = n .÷ 2
    stencil_shape = range.(.-m, m)
    stencil = Array{TIn, dim}(undef, n...)

    function(z_in::AbstractArray{TIn, dim}, z_out::AbstractArray{TOut, dim}, args...)
        @assert (size(z_in) == size(z_out)) "sizes of arrays need to be equal"
        shape = size(z_in)
        for i in CartesianIndices(shape)
            for (k, Δi) in enumerate(CartesianIndices(stencil_shape))
                stencil[k] = offset_value(BT, z_in, i, Δi)
            end
            z_out[i] = f(stencil, args...)
        end
    end
end

stencil(::Type{T}, ::Type{BT}, n::NTuple{dim, Int}, f::Function) where {T, dim, BT <: Boundary{dim}} =
    stencil(T, T, BT, n, f)

convolution(::Type{TIn}, ::Type{TOut}, ::Type{B}, kernel::AbstractArray{U, dim}) where { dim, TIn, TOut, U, B <: Boundary{dim} } =
    stencil(TIn, TOut, B, size(kernel), s -> sum(s .* kernel))

convolution(::Type{B}, kernel::AbstractArray{T, dim}) where {dim, T, B <: Boundary{dim}} =
    stencil(T, T, B, size(kernel), s -> sum(s .* kernel))
# ~/~ end

end
# ~/~ end
