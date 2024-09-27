# ~/~ begin <<docs/src/stencils.md#src/Stencil.jl>>[init]
module Stencil

using ..Config: AbstractBox
using ..BoundaryTrait

export stencil, convolution

# ~/~ begin <<docs/src/stencils.md#stencil-operation>>[init]
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