# ~/~ begin <<docs/src/finite-difference-transport.md#src/Transport/Advection.jl>>[init]
module Advection

using ....CarboKitten: Box
using ...Stencil: stencil!, Size
using ...BoundaryTrait: get_bounded
using ..DifferentialOperators: central_difference, upwind, laplacian
using Unitful
using GeometryBasics
using LinearAlgebra: dot

function advection_coef!(box::Box{BT}, diffusivity, wave_velocity, w, adv, rct) where {BT}
    d = diffusivity
    dx = box.phys_scale
    di = (CartesianIndex(1, 0), CartesianIndex(0, 1))
    for i = eachindex(IndexCartesian(), w)
        v, s = wave_velocity(w[i])

        wx1 = get_bounded(BT, w, i - di[1])
        wx2 = get_bounded(BT, w, i + di[1])
        wy1 = get_bounded(BT, w, i - di[2])
        wy2 = get_bounded(BT, w, i + di[2])

        dw = Vec2((wx2 - wx1) / (2dx), (wy2 - wy1) / (2dx))
        ddw = (wx1 + wx2 + wy1 + wy2 - 4*w[i]) / dx^2

        adv[i] = d * dw + v
        rct[i] = dot(s, dw) - d * ddw
    end
end

function transport_dC!(box::Box{BT}, adv, rct, C, dC) where {BT}
    dx = box.phys_scale
    di = (CartesianIndex(1, 0), CartesianIndex(0, 1))

    @inline upwind(v::T, a, i, di) where {T} =
        if v < zero(T)
            v * (get_bounded(BT, a, i+di) - a[i]) / dx
        else
            v * (a[i] - get_bounded(BT, a, i-di)) / dx
        end

    for i = eachindex(IndexCartesian(), dC)
        dC[i] = rct[i] *  C[i] - upwind(adv[i][1], C, i, di[1]) - upwind(adv[i][2], C, i, di[2]) 
    end

    return dC
end

function max_dt(adv, dx, courant_max)
    u(a) = abs(a[1]) + abs(a[2])
    return courant_max / maximum(u.(adv) ./ dx)
end

"""
    transport!(box, diffusivity, wave_velocity,
               C, w, dC)

Computes `dC` given a `box`, `diffusivity` constant in units of m/Myr,
`wave_velocity` is a function of water depth, returning both velocity in units
of m/Myr, and shear in units of 1/Myr, which should be the derivative of the
velocity w.r.t. water depth. `C` is the concentration of entrained sediment,
`w` the water depth, and `dC` the output derivative of `C`.
"""
function transport!(box::Box{BT}, diffusivity, wave_velocity, C, w, dC) where {BT}
    dx = box.phys_scale
    stencil!(BT, Size(3, 3), dC, C, w) do C, w
        # ~/~ begin <<docs/src/finite-difference-transport.md#advection-transport>>[init]
        d = diffusivity
        v, s = wave_velocity(w[2, 2])

        dw = (central_difference(Val{1}, w, dx), central_difference(Val{2}, w, dx))
        adv = upwind(Val{1}, d * dw[1] + v[1], C, dx) + 
              upwind(Val{2}, d * dw[2] + v[2], C, dx)
        rct = (s[1] * dw[1] + s[2] * dw[2] - d * laplacian(w, dx)) * C[2, 2]
        return rct - adv
        # ~/~ end
    end
end

"""
    transport(box, diffusivity, wave_velocity, wave_shear,
               C, w)

Non-mutating version of [`transport!`](@ref). Allocates and returns `dC`.
"""
function transport(box::Box{BT}, diffusivity, wave_velocity, C::AbstractArray{T}, w) where {BT, T}
	dC = Array{typeof(1.0 * Unitful.unit(T) / u"Myr")}(undef, box.grid_size...)
    transport!(box, diffusivity, wave_velocity, C, w, dC)
    return dC
end

end
# ~/~ end
