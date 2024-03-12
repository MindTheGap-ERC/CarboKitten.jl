# ~/~ begin <<docs/src/transport.md#src/Transport.jl>>[init]
module Transport

using Unitful
using ..Vectors
using ..BoundaryTrait
using ..Config: Box, AbstractBox

# ~/~ begin <<docs/src/boxes.md#vector-offset>>[init]
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
# ~/~ end
# ~/~ begin <<docs/src/transport.md#particle>>[init]
mutable struct Particle{P}
    position::Vec2
    mass::Float64
    critical_stress::Float64
    facies::Int64
    properties::P
end
# ~/~ end
# ~/~ begin <<docs/src/transport.md#particle-transport>>[init]
function transport(box::AbstractBox{BT}, stress, maxit, step) where {BT <: Boundary{2}}
    function (p::Particle{P}) where P
        if p.position ∉ box
            return nothing
        end

        for it in 1:maxit
            if !(p.position ∈ box)
                return nothing
            end
            # @assert (p.position ∈ box) "$(p) in box"
            τ = stress(p)
            if abs(τ) < p.critical_stress
                return p
            end
            Δ = τ * ((box.phys_scale / u"m" |> NoUnits) * step / abs(τ))
            @assert (abs(Δ)*u"m" ≈ box.phys_scale * step) "pos: $(p.position) stress: $(τ) Delta: $(Δ)"
            next_position = offset(box, p.position, Δ)
            if next_position === nothing
                return nothing
            else
                p.position = next_position
            end
        end

        return p
    end
end
# ~/~ end
# ~/~ begin <<docs/src/transport.md#interpolation>>[init]
function interpolate(box::AbstractBox{BT}, f::AbstractMatrix{R}) where {BT <: Boundary{2}, R <: Real}
    p -> interpolate(box, f, p)
end

function interpolate(box::AbstractBox{BT}, f::AbstractMatrix{R}, p::Vec2) where {BT <: Boundary{2}, R <: Real}
    @assert p ∈ box
    phys_scale = box.phys_scale / u"m" |> NoUnits
    l = p / phys_scale
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
         y = (frac.x * (z01 - z00) + (1.0 - frac.x) * (z11 - z10))) / phys_scale

    return (z, ∇)
end

@inline function try_inc(arr, idx, value)
    if idx !== nothing && checkbounds(Bool, arr, idx)
        arr[idx] += value
    end
end

function deposit(box::AbstractBox{BT}, output::Array{Float64,3}) where BT
    function plotter(::Nothing) end

    function plotter(p::Particle{P}) where P
        phys_scale = box.phys_scale / u"m" |> NoUnits
        p.position ∉ box && return
        # we don't want to wrap here. This is done later
        q = p.position - (x=0.5,y=0.5)*phys_scale
        l = q / phys_scale
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
# ~/~ end

end
# ~/~ end