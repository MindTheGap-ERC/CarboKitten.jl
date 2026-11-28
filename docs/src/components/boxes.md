# Box

```component-dag
CarboKitten.Components.Boxes
```

This module makes sure we have access to box properties.

## Box topology

CarboKitten has a 3-dimensional state space, where two dimensions represent cartesian topographic coordinates, and the third dimension is a track record of sedimentation. The cartesian topographic coordinates are always on a regular grid, but depending on the scenario you may choose different map topologies.

- **periodic boundaries** To study sedimentation in a small isolated patch, periodic boundaries seem sufficient. The field is assumed to be infinite in all directions.
- **Von Neumann boundaries** In the case of an island it is nicer to have boundaries with constant derivatives. Produced sediment that flows out of the box is lost to the seas.
- **coastal boundaries** Supposing we simulate a narrow cross section of a coast, we'll have one periodic boundary (in $y$-direction) and Von Neumann boundaries in the $x$-direction.

We parametrize these boundaries as type-level constants in Julia. This way we can use the multiple dispatch mechanism in Julia to obtain specialized implementations for each boundary case, selected at compile time, resulting in efficient run-times.

``` {.julia #boundary-types}
abstract type Boundary{dim} end
struct Reflected{dim} <: Boundary{dim} end
struct Periodic{dim} <: Boundary{dim} end
struct Constant{dim,value} <: Boundary{dim} end
struct Coast <: Boundary{2} end

Base.@deprecate_binding Shelf Coast

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
```

The `Boundary` type is part of the generic `Box` dimension specification.

### Offset indexing

Now we can use these topology types to define three methods for indexing on an offset from some index that is assumed to be within bounds.

``` {.julia #spec}
@testset "offset_value" begin
    @test CartesianIndex(1, 1) == offset_index(Reflected{2}, (3, 3), CartesianIndex(1, 1), CartesianIndex(0, 0))
end
```

``` {.julia #offset-indexing}
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
```

### Canonical coordinates

For both `Periodic` and `Reflected` boundaries it is also possible to write a function that makes any coordinate within bounds. This uses the fact that reflected boundaries are also periodic for a box twice the size.

``` {.julia #canonical-coordinates}
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
```

### Coast boundary

The `Coast` boundary type is specially designed for the simulation of a transect perpendicular to the coast direction. We are periodic in the y-direction and have a Neumannesque constant boundary at the edges of the simulation area.

``` {.julia #offset-indexing}
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
```

### `BoundaryTrait` module

``` {.julia file=src/BoundaryTrait.jl}
# FIXME: Rename this module
module BoundaryTrait

export Boundary, Reflected, Periodic, Constant, Coast, offset_index, offset_value, canonical, get_bounded

<<boundary-types>>
<<offset-indexing>>
<<canonical-coordinates>>

end
```

## Deprecated modules

`Box`, `AbstractBox` and `box_axes` are defined in the `CarboKitten` module itself, and the boundary types live in `BoundaryTrait`. The `Boxes` and `Config` modules used to re-export these names. They are kept only as deprecated aliases, so that existing scripts keep working; import the names from `CarboKitten` instead.

``` {.julia file=src/Boxes.jl}
module Boxes

import ..CarboKitten
import ..BoundaryTrait

Base.@deprecate_binding AbstractBox CarboKitten.AbstractBox
Base.@deprecate_binding Box CarboKitten.Box
Base.@deprecate_binding box_axes CarboKitten.box_axes
Base.@deprecate_binding axes CarboKitten.box_axes false

Base.@deprecate_binding Boundary BoundaryTrait.Boundary false
Base.@deprecate_binding Reflected BoundaryTrait.Reflected false
Base.@deprecate_binding Periodic BoundaryTrait.Periodic false
Base.@deprecate_binding Constant BoundaryTrait.Constant false
Base.@deprecate_binding Coast BoundaryTrait.Coast false
Base.@deprecate_binding Shelf BoundaryTrait.Coast false

end
```

``` {.julia file=src/Config.jl}
module Config

import ..CarboKitten

Base.@deprecate_binding TimeProperties CarboKitten.TimeProperties
Base.@deprecate_binding Box CarboKitten.Box
Base.@deprecate_binding axes CarboKitten.box_axes

end
```

``` {.julia file=test/Components/BoxesSpec.jl}
module BoxesSpec
    using Test
    using CarboKitten.Components.Common
    using CarboKitten.Components.Boxes

    @testset "Components/Boxes" begin
        let box = Box{Periodic{2}}(grid_size=(10, 10), phys_scale=2.0u"m"),
            input = record_type(:Input, Boxes)(box=box)
            @test input.box.grid_size == (10, 10)
            @test input.box.phys_scale == 2.0u"m"
        end
    end
end
```

``` {.julia file=src/Components/Boxes.jl}
@component module Boxes
using ..Common
import ..Common: write_header

@kwdef struct Input <: InputPart
    box::Box
end

function write_header(::Input, input::AbstractInput, output::AbstractOutput)
    x, y = box_axes(input.box)
    set_attribute(output, "x", collect(x) |> in_units_of(u"m"))
    set_attribute(output, "y", collect(y) |> in_units_of(u"m"))
    set_attribute(output, "grid_size_x", input.box.grid_size[1])
    set_attribute(output, "grid_size_y", input.box.grid_size[2])
    set_attribute(output, "phys_scale", input.box.phys_scale |> in_units_of(u"m"))
end
end
```
