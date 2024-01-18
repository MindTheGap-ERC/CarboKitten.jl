# Particle Transport
To compute the transport of material we turn parcels of sediment into particles with coordinates $x, y$, each with their own specific properties $p$, and transport them in a prescribed direction until they reach a point where the shear stress is below critical. Notice that these are not physical particles, rather we choose tracer particles to represent the behaviour of some parcel of sediment.

| property | typical values | units |
|:-------- | --------------:|:----- |
| `grain_size` | ?            | $\mu$m or mm |
| `excess_density` | 1500 | ${\rm kg}/{m^3}$ |
| `critical_stress` | ?     | ${\rm P} = {\rm kg} / ({\rm m} {\rm s}^2)$ |

The `grain_size` and `excess_density` are used in computing the direction a particle will travel. The `critical_stress` determines the stopping criterion.

``` {.julia file=src/Transport.jl}
module Transport

using ..BoundaryTrait

Vec2 = @NamedTuple{x::Float64, y::Float64}
Base.:+(a::Vec2, b::Vec2) = (x=a.x+b.x, y=a.y+b.y)
Base.abs2(a::Vec2) = a.x^2 + a.y^2
Base.abs(a::Vec2) = √(abs2(a))
Base.:*(a::Vec2, b::Float64) = (x=a.x*b, y=a.y*b)
Base.:/(a::Vec2, b::Float64) = (x=a.x/b, y=a.y/b)

struct Particle{P}
    position::Vec2
    critical_stress::Float64
    facies::Int64
    properties::P
end

struct Box
    grid_size::NTuple{2,Int}
    phys_size::Vec2
    phys_scale::Float64
end

function offset(::Type{Reflected{2}}, box::Box, a::Vec2, Δa::Vec2)
    clip(i, a, b) = (i < a ? a + a - i : (i > b ? b + b - i : i))
    (x=clip(a.x+Δa.x, 0.0, box.phys_size.x)
    ,y=clip(a.y+Δa.y, 0.0, box.phys_size.y))
end

function offset(::Type{Periodic{2}}, box::Box, a::Vec2, Δa::Vec2)
    (x=(a.x+Δa.x) % box.phys_size.x
    ,y=(a.y+Δa.y) % box.phys_size.y)
end

function offset(::Type{Constant{2,Value}}, box::Box, a::Vec2, Δa::Vec2)
    b = a + Δa
    if b.x < 0.0 | b.x > box.phys_size.x | b.y < 0.0 | b.y > box.phys_size.y
        nothing
    else
        b
    end
end

function transport(::Type{BT}, box::Box, Δx::Float64, stress, p::Particle{P}) where {BT <: Boundary{2}, P}
    while true
        τ = stress(p)
        if abs(τ) > p.critical_stress
            return p
        end
        Δ = τ * (Δx / abs(τ))
        p.position = offset(BT, box, p.position, Δ)
        if p.position === nothing
            return nothing
        end
    end
end

function deposit()
end

end
```
