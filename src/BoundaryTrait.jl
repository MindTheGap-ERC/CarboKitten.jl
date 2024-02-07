# ~/~ begin <<docs/src/stencils.md#src/BoundaryTrait.jl>>[init]
module BoundaryTrait

export Boundary, Reflected, Periodic, Constant, Shelf, offset_index, offset_value

# ~/~ begin <<docs/src/stencils.md#boundary-types>>[init]
abstract type Boundary{dim} end
struct Reflected{dim} <: Boundary{dim} end
struct Periodic{dim} <: Boundary{dim} end
struct Constant{dim,value} <: Boundary{dim} end
# ~/~ end
# ~/~ begin <<docs/src/stencils.md#boundary-types>>[1]
struct Shelf <: Boundary{2} end
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
# ~/~ begin <<docs/src/stencils.md#offset-indexing>>[1]
function offset_index(::Type{Shelf}, shape::NTuple{2,Int}, i::CartesianIndex, Δi::CartesianIndex)
    j = i + Δi
    j[1] >= 1 | j[1] <= shape[1] ? CartesianIndex(j[1], mod1(j[2], shape[2])) : nothing
end

function offset_value(::Type{Shelf}, z::AbstractArray, i::CartesianIndex, Δi::CartesianIndex)
    j = i + Δi
    if j[1] < 1
        return z[1, mod1(j[2], shape[2])]
    elseif j[1] > shape[1]
        return z[shape[1], mod1(j[2], shape[2])]
    else
        return z[j[1], mod1(j[2], shape[2])]
    end
end
# ~/~ end
# ~/~ begin <<docs/src/stencils.md#canonical-coordinates>>[init]
function canonical(::Type{Periodic{2}}, shape::NTuple{dim,Int}, i::CartesianIndex) where {dim}
    CartesianIndex(mod1.(Tuple(i), shape)...)
end

function canonical(::Type{Reflected{2}}, shape::NTuple{dim,Int}, i::CartesianIndex) where {dim}
    modflip(a, l) = let b = mod1(a, 2l)
        b > l ? 2l - b : b
    end
    CartesianIndex(modflip.(Tuple(i), shape)...) 
end
# ~/~ end

end
# ~/~ end