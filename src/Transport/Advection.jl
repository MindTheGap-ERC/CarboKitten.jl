# ~/~ begin <<docs/src/finite-difference-transport.md#src/Transport/Advection.jl>>[init]
module Advection

using ....CarboKitten: Box
using ...Stencil: stencil!, Size
using ..DifferentialOperators: central_difference, upwind, laplacian
using Unitful

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
function transport(box::Box{BT}, diffusivity, wave_velocity, C, w) where {BT}
    dC = Array{typeof(1.0/u"yr")}(undef, box.grid_size...)
    transport!(box, diffusivity, wave_velocity, C, w, dC)
    return dC
end

end
# ~/~ end
