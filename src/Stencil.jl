# ~/~ begin <<docs/src/stencils.md#src/Stencil.jl>>[init]
module Stencil

export Boundary, Reflected, Periodic, Constant, stencil, convolution, offset_index, offset_value

# ~/~ begin <<docs/src/stencils.md#boundary-trait>>[init]
abstract type Boundary{dim} end
struct Reflected{dim} <: Boundary{dim} end
struct Periodic{dim} <: Boundary{dim} end
struct Constant{dim,value} <: Boundary{dim} end
# ~/~ end
# ~/~ begin <<docs/src/stencils.md#offset-indexing>>[init]
function offset_index(::Type{Periodic{dim}}, shape::NTuple{dim,Int}, i::CartesianIndex, Δi::CartesianIndex) where {dim}
    CartesianIndex(mod1.(Tuple(i + Δi), shape)...)
end

function offset_index(::Type{Reflected{dim}}, shape::NTuple{dim,Int}, i::CartesianIndex, Δi::CartesianIndex) where {dim}
    clip(i, a, b) = (i < a ? a + a - i : (i > b ? b + b - i : i))
    CartesianIndex(clip.(Tuple(i + Δi), ones(Int, dim), shape)...)
end

function offset_index(::Type{Constant{dim,value}}, shape::NTuple{dim,Int}, i::CartesianIndex, Δi::CartesianIndex) where {dim,value}
    j = i + Δi
    all(checkindex.(Bool, range.(1, shape), Tuple(j))) ? j : nothing
end

function offset_value(BT::Type{B}, z::AbstractArray, i::CartesianIndex, Δi::CartesianIndex) where {dim,B<:Boundary{dim}}
    z[offset_index(BT, size(z), i, Δi)]
end

function offset_value(::Type{Constant{dim,value}}, z::AbstractArray, i::CartesianIndex, Δi::CartesianIndex) where {dim,value}
    j = i + Δi
    (checkbounds(Bool, z, j) ? z[j] : value)
end
# ~/~ end
# ~/~ begin <<docs/src/stencils.md#stencil-operation>>[init]
function stencil(::Type{T}, ::Type{BT}, n::NTuple{dim,Int}, f::Function) where {T,dim,BT<:Boundary{dim}}
    return stencil(T, T, BT, n, f)
end

function stencil(::Type{In}, ::Type{Out}, ::Type{BT}, n::NTuple{dim,Int}, f::Function) where {In,Out,dim,BT<:Boundary{dim}}
    m = n .÷ 2
    stencil_shape = range.(.-m, m)
    stencil = zeros(T, n)

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