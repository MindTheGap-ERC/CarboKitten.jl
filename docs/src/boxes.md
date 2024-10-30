# Generic Parameters

``` {.julia file=src/Config.jl}
module Config

export TimeProperties

using ..BoundaryTrait
using ..Vectors

using Unitful
using Unitful.DefaultSymbols

<<config-types>>

end
```

Physical parameters of CarboKitten all should have units, see our [refresher on `Unitful.jl`](unitful.md).

## Time properties

Time stepping is specified in `TimeProperties`. We'll have `time_steps` number of time steps, each of physical time `Δt`. However, only one in `write_interval` steps is written to disk.

``` {.julia #config-types}
abstract type AbstractTimeProperties end

@kwdef struct TimeProperties <: AbstractTimeProperties
    t0::typeof(1.0u"Myr") = 0.0u"Myr"
    Δt::typeof(1.0u"Myr")
    steps::Int
    write_interval::Int = 1
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
