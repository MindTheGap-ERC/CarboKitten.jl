# ~/~ begin <<docs/src/components/boxes.md#src/BoundaryTrait.jl>>[init]
# FIXME: Rename this module
module BoundaryTrait

export Boundary, Reflected, Periodic, Constant, Coast, Shelf, offset_index, offset_value, canonical, get_bounded

# ~/~ begin <<docs/src/components/boxes.md#boundary-types>>[init]
abstract type Boundary{dim} end
struct Reflected{dim} <: Boundary{dim} end
struct Periodic{dim} <: Boundary{dim} end
struct Constant{dim,value} <: Boundary{dim} end
struct Coast <: Boundary{2} end

const Shelf = Coast  # FIXME: Old name, should be removed

@inline modflip(a, l) = let b = mod1(a, 2l)
    b > l ? 2l - b + 1 : b
end

@inline get_bounded(::Type{Constant{dim, value}}, a, i) where {dim, value} =
    checkbounds(Bool, a, i) ? a[i] : val

@inline get_bounded(::Type{Periodic{dim}}, a, i) where {dim} =
    checkbounds(Bool, a, i) ? a[i] : a[mod1.(Tuple(i), size(a))...]

@inline get_bounded(::Type{Reflected{dim}}, a, i) where {dim} =
    checkbounds(Bool, a, i) ? a[i] : a[modflip.(Tuple(i), size(a))...]

@inline get_bounded(::Type{Coast}, a, i) =
    checkbounds(Bool, a, i) ? a[i] : a[modflip(i[1], size(a)[1]), mod1(i[2], size(a)[2])]
# ~/~ end
# ~/~ begin <<docs/src/components/boxes.md#offset-indexing>>[init]
function offset_index(::Type{BT}, shape::NTuple{dim,Int}, i::CartesianIndex, Δi::CartesianIndex) where {dim, BT <: Boundary{dim}}
    canonical(BT, shape, i + Δi)
end

function offset_value(BT::Type{B}, z::AbstractArray, i::CartesianIndex, Δi::CartesianIndex) where {dim,B<:Boundary{dim}}
    z[offset_index(BT, size(z), i, Δi)]
end

function offset_value(::Type{Constant{dim,value}}, z::AbstractArray, i::CartesianIndex, Δi::CartesianIndex) where {dim,value}
    j = i + Δi
    (checkbounds(Bool, z, j) ? z[j] : value)
end
# ~/~ end
# ~/~ begin <<docs/src/components/boxes.md#offset-indexing>>[1]
function canonical(::Type{Coast}, shape::NTuple{2, Int}, i::CartesianIndex)
    if i[1] < 1 || i[1] > shape[1]
        return nothing
    end
    return CartesianIndex(i[1], mod1(i[2], shape[2]))
end

function offset_value(::Type{Coast}, z::AbstractArray, i::CartesianIndex, Δi::CartesianIndex)
    j = i + Δi
    shape = size(z)
    if j[1] < 1
        return z[1, mod1(j[2], shape[2])]
    elseif j[1] > shape[1]
        return z[shape[1], mod1(j[2], shape[2])]
    else
        return z[j[1], mod1(j[2], shape[2])]
    end
end
# ~/~ end
# ~/~ begin <<docs/src/components/boxes.md#canonical-coordinates>>[init]
function canonical(::Type{Periodic{dim}}, shape::NTuple{dim,Int}, i::CartesianIndex) where {dim}
    CartesianIndex(mod1.(Tuple(i), shape)...)
end

function canonical(::Type{Reflected{dim}}, shape::NTuple{dim,Int}, i::CartesianIndex) where {dim}
    modflip(a, l) = let b = mod1(a, 2l)
        b > l ? 2l - b + 1 : b
    end
    CartesianIndex(modflip.(Tuple(i), shape)...) 
end

function canonical(::Type{Constant{dim, value}}, shape::NTuple{dim,Int}, i::CartesianIndex) where {dim, value}
    all(checkindex.(Bool, range.(1, shape), Tuple(i))) ? i : nothing
end
# ~/~ end

end
# ~/~ end