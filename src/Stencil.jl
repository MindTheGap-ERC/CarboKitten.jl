# ~/~ begin <<docs/src/stencils.md#src/Stencil.jl>>[init]
module Stencil

using ..BoundaryTrait
export stencil, convolution

# ~/~ begin <<docs/src/stencils.md#stencil-operation>>[init]
function stencil(::Type{T}, ::Type{BT}, n::NTuple{dim,Int}, f::Function) where {T,dim,BT<:Boundary{dim}}
    return stencil(T, T, BT, n, f)
end

function stencil(::Type{In}, ::Type{Out}, ::Type{BT}, n::NTuple{dim,Int}, f::Function) where {In,Out,dim,BT<:Boundary{dim}}
    m = n .÷ 2
    stencil_shape = range.(.-m, m)
    stencil = zeros(In, n)

    function (z_in::AbstractArray{In,dim}, z_out::AbstractArray{Out,dim}, args...)
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

function convolution(::Type{B}, kernel::Array{T,dim}) where {dim,T,B<:Boundary{dim}}
    stencil(T, B, size(kernel), s -> sum(s .* kernel))
end
# ~/~ end

end
# ~/~ end