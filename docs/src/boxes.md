# Generic Parameters

``` {.julia file=src/Config.jl}
module Config

export AbstractBox, Box

using ..BoundaryTrait
using ..Vectors

using Unitful
using Unitful.DefaultSymbols

<<config-types>>

end
```

Physical parameters of CarboKitten all should have units, see our [refresher on `Unitful.jl`](unitful.md).

## Box topology

CarboKitten has a 3-dimensional state space, where two dimensions represent cartesian topographic coordinates, and the third dimension is a track record of sedimentation. The cartesian topographic coordinates are always on a regular grid, but depending on the scenario you may choose different map topologies.

- **periodic boundaries** To study sedimentation in a small isolated patch, periodic boundaries seem sufficient. The field is assumed to be infinite in all directions.
- **Von Neumann boundaries** In the case of an island it is nicer to have boundaries with constant derivatives. Produced sediment that flows out of the box is lost to the seas.
- **coastal shelf boundaries** Supposing we simulate a narrow cross section of a carbonate shelf, we'll have one periodic boundary (in $y$-direction) and Von Neumann boundaries in the $x$-direction.

We parametrize these boundaries as type-level constants in Julia. This way we can use the multiple dispatch mechanism in Julia to obtain specialized implementations for each boundary case, selected at compile time, resulting in efficient run-times.

``` {.julia #boundary-types}
abstract type Boundary{dim} end
struct Reflected{dim} <: Boundary{dim} end
struct Periodic{dim} <: Boundary{dim} end
struct Constant{dim,value} <: Boundary{dim} end
struct Shelf <: Boundary{2} end
```

The `Boundary` type is part of the generic `Box` dimension specification.

``` {.julia #config-types}
abstract type AbstractBox{BT} end

struct Box{BT} <: AbstractBox{BT}
    grid_size::NTuple{2,Int}
    phys_scale::typeof(1.0m)
    phys_size::Vec2

    function Box{BT}(;grid_size::NTuple{2, Int}, phys_scale::Quantity{Float64, 𝐋, U}) where {BT <: Boundary{2}, U}
        new{BT}(grid_size, phys_scale, phys_size(grid_size, phys_scale))
    end
end

function axes(box::Box)
	y_axis = (0:(box.grid_size[2] - 1)) .* box.phys_scale
	x_axis = (0:(box.grid_size[1] - 1)) .* box.phys_scale
	return x_axis, y_axis
end

phys_size(grid_size, phys_scale) = (
    x = grid_size[1] * (phys_scale / m |> NoUnits),
    y = grid_size[2] * (phys_scale / m |> NoUnits))
```

Now we can specify the box parameters as follows:

``` {.julia file=test/ConfigSpec.jl}
@testset "Config" begin
    using CarboKitten.BoundaryTrait
    using CarboKitten.Config: Box
    using CarboKitten.Vectors

    box = Box{Shelf}(
        grid_size = (100, 50),
        phys_scale = 1.0u"km")
    @test box.phys_size == (x=100000.0, y=50000.0)
end
```

## Time properties

Time stepping is specified in `TimeProperties`. We'll have `time_steps` number of time steps, each of physical time `Δt`. However, only one in `write_interval` steps is written to disk.

``` {.julia #config-types}
abstract type AbstractTimeProperties end

@kwdef struct TimeProperties <: AbstractTimeProperties
    Δt::typeof(1.0u"Myr")
    steps::Int
    write_interval::Int
end
```

## Vectors

To trace the position of particles we define a `NamedTuple` with `x` and `y` members and define common vector operations on those.

``` {.julia file=src/Vectors.jl}
module Vectors

export Vec2

Vec2 = @NamedTuple{x::Float64, y::Float64}
Base.:+(a::Vec2, b::Vec2) = (x=a.x+b.x, y=a.y+b.y)
Base.abs2(a::Vec2) = a.x^2 + a.y^2
Base.abs(a::Vec2) = √(abs2(a))
Base.:*(a::Vec2, b::Float64) = (x=a.x*b, y=a.y*b)
Base.:/(a::Vec2, b::Float64) = (x=a.x/b, y=a.y/b)
Base.:*(a::Float64, b::Vec2) = b*a
Base.:-(a::Vec2, b::Vec2) = (x=a.x-b.x, y=a.y-b.y)
Base.:-(a::Vec2) = (x=-a.x, y=-a.y)
Base.zero(::Type{Vec2}) = (x=0.0, y=0.0)

end
```

## Offset indexing

Now we can use these traits to define three methods for indexing on an offset from some index that is assumed to be within bounds.

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

``` {.julia file=src/BoundaryTrait.jl}
module BoundaryTrait

export Boundary, Reflected, Periodic, Constant, Shelf, offset_index, offset_value, canonical

<<boundary-types>>
<<offset-indexing>>
<<canonical-coordinates>>

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

### Shelf boundary

The `Shelf` boundary type is specially designed for the simulation of a transect perpendicular to the coast direction. We are periodic in the y-direction and have a Neumannesque constant boundary at the edges of the simulation area.

``` {.julia #offset-indexing}
function canonical(::Type{Shelf}, shape::NTuple{2, Int}, i::CartesianIndex)
    if i[1] < 1 || i[1] > shape[1]
        return nothing
    end
    return CartesianIndex(i[1], mod1(i[2], shape[2]))
end

function offset_value(::Type{Shelf}, z::AbstractArray, i::CartesianIndex, Δi::CartesianIndex)
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

## Boxes

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
