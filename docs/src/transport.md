# Particle Transport
To compute the transport of material we turn parcels of sediment into particles with coordinates $x, y$, each with their own specific properties $p$, and transport them in a prescribed direction until they reach a point where the shear stress is below critical. Notice that these are not physical particles, rather we choose tracer particles to represent the behaviour of some parcel of sediment.

| property | typical values | units |
|:-------- | --------------:|:----- |
| `grain_size` | ?            | $\mu$m or mm |
| `excess_density` | 1500 | ${\rm kg}/{m^3}$ |
| `critical_stress` | ?     | ${\rm P} = {\rm kg} / ({\rm m} {\rm s}^2)$ |

The `grain_size` and `excess_density` are used in computing the direction a particle will travel. The `critical_stress` determines the stopping criterion.

## Particles
In this section we only deal with those properties of particles that are required to model their transport and (re)deposition. We leave the computation of the stress for later. 

!!! note "Abstraction and Inheritance"

    We don't know all the properties yet that might be involved in computing that stress, so we'll make the `Particle` type extensible. The classic way to do so in computer science is inheritance, but Julia doesn't support inheritance. Instead, you might define an `abstract type` and define methods for that type, then require derived `struct`s to have the same properties, or overload the methods. However, that is not so nice and actually error prone. Here I chose to make the `Particle` type extensible by using a type argument for the `properties` member. In the simplest case you can define an alias `const Particle = Transport.Particle{Nothing}`.

We need the position of a particle to track its orbit, its mass and facies type to compute the deposition onto a grid and the critical stress to tell when the particle should stop moving.

``` {.julia #particle}
mutable struct Particle{P}
    position::Vec2
    mass::Float64
    critical_stress::Float64
    facies::Int64
    properties::P
end
```

The choice is to have the `position` member represent coordinates in logical (pixel based) units. An alternative would be to use logical (pixel based) coordinates, possibly saving a few floating-point operations. The hope is that using physical coordinates saves some mental capacity in remembering whether coordinates were converted or not, but this is not a strong preference. Do notice that a user wanting to tinker with transport is exposed to this API.

## Interpolation and Deposition
Bi-linear interpolation and mass deposition are in some ways very similar problems. We deposit a particle as a little square shaped object on the grid. We compute the area of the intersection of the pixel for a group of $2\times 2$ pixels. This is the same weighting scheme used in bi-linear interpolation. This particular interpolation routine also returns the local gradient.

### Specs
Deposition should be mass conserving, so we test for a number of particles that the total adds up.

``` {.julia #transport-spec}
using CarboKitten.Transport: deposit, Box, Particle
using CarboKitten.BoundaryTrait: Periodic

TestParticle = Particle{Nothing}

box = Box((16, 16), (x=1.0, y=1.0), 1.0/16.0)
n_particles = 42
particles = rand(Float64, 2, n_particles) |> eachcol .|> 
		(p -> TestParticle((x=p[1], y=p[2]), 1.0, 1.0, 1, nothing));

density = let
	target = zeros(Float64, 1, 16, 16)
	particles .|> deposit(Periodic{2}, box, target)
	target
end
@test sum(density) ≈ n_particles
```

Then we test two special cases: one where the particle is exactly at position $0, 0$, and the other at $0.5, 0.5$ (pixel coordinates).

``` {.julia #transport-spec}
density = let
    target = zeros(Float64, 1, 16, 16)
    Particle((x=0.0, y=0.0), 1.0, 1.0, 1, nothing) |> deposit(Periodic{2}, box, target)
    target
end
@test density[1,1,1] ≈ 0.25
@test density[1,1,end] ≈ 0.25
@test density[1,end,1] ≈ 0.25
@test density[1,end,end] ≈ 0.25
```

``` {.julia #transport-spec}
density = let
    target = zeros(Float64, 1, 16, 16)
    Particle((x=box.phys_scale/2, y=box.phys_scale/2),
        1.0, 1.0, 1, nothing) |> deposit(Periodic{2}, box, target)
    target
end
@test density[1,1,1] ≈ 1.0
@test density[1,1,end] ≈ 0
@test density[1,end,1] ≈ 0
@test density[1,end,end] ≈ 0
```

``` {.julia #transport-spec}
height = [0.0 0.0; 1.0 1.0]
```

### Implementation

``` {.julia #interpolation}
function interpolate(::Type{BT}, box::Box, f::AbstractMatrix{R}) where {BT <: Boundary{2}, R <: Real}
    p -> interpolate(BT, box, f, p)
end

function interpolate(::Type{BT}, box::Box, f::AbstractMatrix{R}, p::Vec2) where {BT <: Boundary{2}, R <: Real}
    @assert p ∈ box
    l = p / box.phys_scale
    node = (x=floor(l.x) + 1.0, y=floor(l.y) + 1.0)
    idx = CartesianIndex(Int(node.x), Int(node.y))
    frac = node - l
    z00 = f[idx]
    z10 = offset_value(BT, f, idx, CartesianIndex(1, 0))
    z01 = offset_value(BT, f, idx, CartesianIndex(0, 1))
    z11 = offset_value(BT, f, idx, CartesianIndex(1, 1))
    z = frac.x * frac.y * z00 + (1.0 - frac.x) * frac.y * z10 +
        frac.x * (1.0 - frac.y) * z01 + (1.0 - frac.x) * (1.0 - frac.y) * z11
    ∇ = (x = (frac.y * (z10 - z00) + (1.0 - frac.y) * (z11 - z01)),
         y = (frac.x * (z01 - z00) + (1.0 - frac.x) * (z11 - z10))) / box.phys_scale

    return (z, ∇)
end

@inline function try_inc(arr, idx, value)
    if idx !== nothing && checkbounds(Bool, arr, idx)
        arr[idx] += value
    end
end

function deposit(::Type{BT}, box::Box, output::Array{Float64,3}) where BT
    function plotter(::Nothing) end

    function plotter(p::Particle{P}) where P
        p.position ∉ box && return
        # we don't want to wrap here. This is done later
        q = p.position - (x=-0.5,y=-0.5)*box.phys_scale
        l = q / box.phys_scale
        node = (x=floor(l.x) + 1.0, y=floor(l.y) + 1.0)
        idx = CartesianIndex(Int(node.x), Int(node.y))
        frac = node - l
        slice = @view output[p.facies,:,:]
        try_inc(slice, canonical(BT, box.grid_size, idx), frac.x * frac.y * p.mass)
        try_inc(slice, offset_index(BT, box.grid_size, idx, CartesianIndex(0, 1)),
            frac.x * (1.0 - frac.y) * p.mass)
        try_inc(slice, offset_index(BT, box.grid_size, idx, CartesianIndex(1, 0)),
            (1.0 - frac.x) * frac.y * p.mass)
        try_inc(slice, offset_index(BT, box.grid_size, idx, CartesianIndex(1, 1)),
            (1.0 - frac.x) * (1.0 - frac.y) * p.mass)
    end

    plotter
end
```

## Transport

``` {.julia #particle-transport}
function transport(::Type{BT}, box::Box, stress, maxit, step) where {BT <: Boundary{2}}
    function (p::Particle{P}) where P
        if p.position ∉ box
            return nothing
        end

        for it in 1:maxit
            @assert (p.position ∈ box) "$(p) in box"
            τ = stress(p)
            if abs(τ) < p.critical_stress
                return p
            end
            Δ = τ * (box.phys_scale * step / abs(τ))
            @assert (abs(Δ) ≈ box.phys_scale * step) "pos: $(p.position) stress: $(τ) Delta: $(Δ)"
            next_position = offset(BT, box, p.position, Δ)
            if next_position === nothing
                return nothing
            else
                p.position = next_position
            end
        end

        return p
    end
end
```



``` {.julia file=src/Transport.jl}
module Transport

using ..Vectors
using ..BoundaryTrait

<<vector-offset>>
<<particle>>
<<particle-transport>>
<<interpolation>>

end
```

``` {.julia file=test/TransportSpec.jl}
@testset "TransportSpec" begin
    <<transport-spec>>
end
```
