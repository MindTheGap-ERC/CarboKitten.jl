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

const Shelf = Coast  # FIXME: Old name, should be removed
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

export Boundary, Reflected, Periodic, Constant, Coast, Shelf, offset_index, offset_value, canonical

<<boundary-types>>
<<offset-indexing>>
<<canonical-coordinates>>

end
```

## Box properties

``` {.julia #box-type}
const axes = box_axes

phys_size(grid_size, phys_scale) = (
    x = grid_size[1] * (phys_scale / m |> NoUnits),
    y = grid_size[2] * (phys_scale / m |> NoUnits))
```

### Floating-point vectors

We need to define how particles move past boundaries. Similar to the grid based `offset_index` method, we define the `offset` method for a `Vec2`.

``` {.julia #vector-offset}
Base.in(a::Vec2, box::Box) =
    a.x >= 0.0 && a.x < box.phys_size.x && a.y >= 0.0 && a.y < box.phys_size.y

function offset(box::AbstractBox{Reflected{2}}, a::Vec2, Δa::Vec2)
    clip(i, a, b) = (i < a ? a + a - i : (i > b ? b + b - i : i))
    (x=clip(a.x+Δa.x, 0.0, box.phys_size.x)
    ,y=clip(a.y+Δa.y, 0.0, box.phys_size.y))
end

function offset(box::AbstractBox{Periodic{2}}, a::Vec2, Δa::Vec2)
    (x=mod(a.x+Δa.x, box.phys_size.x)
    ,y=mod(a.y+Δa.y, box.phys_size.y))
end

function offset(box::AbstractBox{Constant{2,Value}}, a::Vec2, Δa::Vec2) where Value
    b = a + Δa
    if b ∉ box
        nothing
    else
        b
    end
end

function offset(box::AbstractBox{Shelf}, a::Vec2, Δa::Vec2)
    b = a + Δa
    if b.x < 0.0 || b.x >= box.phys_size.x
        nothing
    else
        (x=b.x, y=mod(b.y, box.phys_size.y))
    end
end
```

## Modules

``` {.julia file=src/Boxes.jl}
module Boxes

using ..BoundaryTrait
using ..CarboKitten: AbstractBox, Box, box_axes
using ..Vectors

using Unitful
using Unitful.DefaultSymbols

export AbstractBox, Box, box_axes

<<box-type>>
<<vector-offset>>

end
```

``` {.julia file=test/Components/BoxesSpec.jl}
module BoxesSpec
    using Test
    using CarboKitten.Components.Common
    using CarboKitten.Components.Boxes

    @testset "Components/Boxes" begin
        let box = Box{Periodic{2}}(grid_size=(10, 10), phys_scale=2.0u"m"),
            input = Boxes.Input(box=box)
            @test input.box.grid_size == (10, 10)
            @test input.box.phys_scale == 2.0u"m"
        end
    end
end
```

``` {.julia file=src/Components/Boxes.jl}
@compose module Boxes
using ..Common

@kwdef struct Input <: AbstractInput
    box::Box
end

function write_header(fid, input::AbstractInput)
    x, y = box_axes(input.box)

    gid = fid["input"]
    gid["x"] = collect(x) |> in_units_of(u"m")
    gid["y"] = collect(y) |> in_units_of(u"m")
end
end
```
