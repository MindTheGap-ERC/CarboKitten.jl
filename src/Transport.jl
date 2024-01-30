# ~/~ begin <<docs/src/transport.md#src/Transport.jl>>[init]
module Transport

using ..Vectors
using ..BoundaryTrait

struct Particle{P}
    position::Vec2
    mass::Float64
    critical_stress::Float64
    facies::Int64
    properties::P
end

struct Box
    grid_size::NTuple{2,Int}
    phys_size::Vec2
    phys_scale::Float64
    n_facies::Int
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

function offset(::Type{Constant{2,Value}}, box::Box, a::Vec2, Δa::Vec2) where Value
    b = a + Δa
    if b.x < 0.0 | b.x > box.phys_size.x | b.y < 0.0 | b.y > box.phys_size.y
        nothing
    else
        b
    end
end

function offset(::Type{Shelf}, box::Box, a::Vec2, Δa::Vec2)
    b = a + Δa
    if b.x < 0.0 | b.x > box.phys_size.x
        return nothing
    else
        b.y = b.y % box.phys_size.y
        return b
    end
end

function interpolate(::Type{BT}, box::Box, f::AbstractMatrix{R}, p::Vec2) where R <: Real
    node = (x=ceil(p.x / box.phys_scale), y=ceil(p.y / box.phys_scale))
    idx = CartesianIndex(Int(node.x), Int(node.y))
    frac = (x=node.x - p.x / box.phys_scale, y=node.y - p.y / box.phys_scale)
    z00 = f[idx]
    z10 = offset_value(BT, f, idx, CartesianIndex(1, 0))
    z01 = offset_value(BT, f, idx, CartesianIndex(0, 1))
    z11 = offset_value(BT, f, idx, CartesianIndex(1, 1))
    z = frac.x * frac.y * z00 + (1.0 - frac.x) * frac.y * z10 +
        frac.x * (1.0 - frac.y) * z01 + (1.0 - frac.x) * (1.0 - frac.y) * z11
    ∇ = (x = (frac.y * (z10 - z00) + (1.0 - frac.y) * (z11 - z01)) / box.phys_scale,
         y = (frac.x * (z01 - z00) + (1.0 - frac.x) * (z11 - z10)) / box.phys_scale)

    return (z, ∇)
end

function transport(::Type{BT}, box::Box, stress) where {BT <: Boundary{2}}
    function (p::Particle{P}) where P
        while true
            τ = stress(p)
            if abs(τ) > p.critical_stress
                return p
            end
            Δ = τ * (box.phys_scale / abs(τ))
            next_position = offset(BT, box, p.position, Δ)
            if next_position === nothing
                return nothing
            else
                p.position = next_position
            end
        end
    end
end

function deposit(::Type{BT}, box::Box, output::Array{Float64,3}) where BT
    function (p::Particle{P}) where P
        node = (x=ceil(p.x / box.phys_scale), y=ceil(p.y / box.phys_scale))
        idx = CartesianIndex(Int(node.x), Int(node.y))
        frac = (x=node.x - p.x / box.phys_scale, y=node.y - p.y / box.phys_scale)
        slice = @view output[p.facies,:,:]
        slice[idx] += frac.x * frac.y * p.mass
        slice[offset_index(BT, box.grid_size, idx, CartesianIndex(0, 1))] +=
            frac.x * (1.0 - frac.y) * p.mass
        slice[offset_index(BT, box.grid_size, idx, CartesianIndex(1, 0))] +=
            (1.0 - frac.x) * frac.y * p.mass
        slice[offset_index(BT, box.grid_size, idx, CartesianIndex(1, 1))] +=
            (1.0 - frac.x) * (1.0 - frac.y) * p.mass
    end
end

end
# ~/~ end