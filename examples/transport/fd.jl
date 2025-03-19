# ~/~ begin <<docs/src/finite-difference-transport.md#examples/transport/fd.jl>>[init]
using ModuleMixins
using Unitful
using CarboKitten

module Runner
using CarboKitten
using ProgressLogging

n_steps(input) = input.time.steps

function run_model(f, ::Type{Model{M}}, input) where {M}
    state = M.initial_state(input)
    f(0, state)

    @progress for i = 1:n_steps(input)
        M.step!(input, state)
        f(i, state)
    end

    return state
end

do_nothing(_i, _s) = ()

run_model(::Type{Model{M}}, input) where {M} = run_model(do_nothing, Model{M}, input)
end


module DifferentialOperators

laplacian(a::AbstractMatrix, dx) =
    (a[1, 2] + a[2, 1] + a[3, 2] + a[2, 3] - 4 * a[2, 2]) / dx^2

central_difference(::Type{Val{1}}, a::AbstractMatrix, dx) =
    (a[3, 2] - a[1, 2]) / (2dx)

central_difference(::Type{Val{2}}, a::AbstractMatrix, dx) =
    (a[2, 3] - a[2, 1]) / (2dx)

upwind(v::T, a1, a2, a3, dx) where {T} =
    if v < zero(T)
        v * (a3 - a2) / dx
    else
        v * (a2 - a1) / dx
    end

upwind(::Type{Val{1}}, v, a::AbstractMatrix, dx) =
    upwind(v, a[:, 2]..., dx)

upwind(::Type{Val{2}}, v, a::AbstractMatrix, dx) =
    upwind(v, a[2, :]..., dx)

end

module Solvers
using Unitful

function runge_kutta_4(box)
    k1 = Array{typeof(1.0u"1/yr")}(undef, box.grid_size...)
    k2 = Array{typeof(1.0u"1/yr")}(undef, box.grid_size...)
    k3 = Array{typeof(1.0u"1/yr")}(undef, box.grid_size...)
    k4 = Array{typeof(1.0u"1/yr")}(undef, box.grid_size...)
    function (df, y, t, dt)
        k1 .= df(y, t)
        k2 .= df(y .+ dt/2 .* k1, t + dt/2)
        k3 .= df(y .+ dt/2 .* k2, t + dt/2)
        k4 .= df(y .+ dt .* k3, t + dt)
        y .+= (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4) .* (dt/6)
    end
end

function forward_euler(df, y, t, dt)
    y .+= dt .* df(y, t)
end

end

@compose module Advection

using CarboKitten  #: Box, TimeProperties
using CarboKitten.Stencil: stencil!, Size
using Unitful
using ..DifferentialOperators: upwind

@kwdef struct Input
    box::Box
    time::TimeProperties
    initial_state::Array{Float64}

    wave_velocity
    solver
end

@kwdef mutable struct State
    time::typeof(1.0u"Myr")
    # delta::Array{Float64}
    value::Array{Float64}
end

initial_state(input) = State(
    time = input.time.t0,
    # delta = zeros(Float64, input.box.grid_size...),
    value = copy(input.initial_state))

advection(v, a::AbstractMatrix, dx) =
    - (upwind(Val{1}, v[1], a, dx) + upwind(Val{2}, v[2], a, dx))

pure_advection!(box::Box{BT}, v, a, da) where {BT} = stencil!(BT, Size(3, 3), da, a) do a
    advection(v, a, box.phys_scale)
end

pure_advection(box::Box{BT}, v, a) where {BT} = let da = Array{typeof(1.0/u"yr")}(undef, box.grid_size...)
    pure_advection!(box, v, a, da)
    return da
end

function step!(input, state)
    input.solver(
        (a, t) -> pure_advection(input.box, input.wave_velocity, a),
        state.value, state.time, input.time.Δt)
    state.time += input.time.Δt
end

end

module Diffusion

using CarboKitten  #: Box, TimeProperties
using CarboKitten.Stencil: stencil!, Size
using Unitful
using ..DifferentialOperators: central_difference, upwind, laplacian

@kwdef struct Input
    box::Box
    time::TimeProperties
    initial_state::Array{Float64}
    topography::Array{typeof(1.0u"m")}

    diffusion_coefficient
    solver
end

@kwdef mutable struct State
    time::typeof(1.0u"Myr")
    # delta::Array{Float64}
    value::Array{Float64}
end

initial_state(input) = State(
    time = input.time.t0,
    # delta = zeros(Float64, input.box.grid_size...),
    value = copy(input.initial_state))

function diffusion!(box::Box{BT}, d, C, w, dC) where {BT}
    dx = box.phys_scale
    stencil!(BT, Size(3, 3), dC, C, w) do C, w
        dw = (central_difference(Val{1}, w, dx), central_difference(Val{2}, w, dx))
        diff = -d * C[2, 2] * laplacian(w, dx)
        adve = -d * (upwind(Val{1}, dw[1], C, dx) + upwind(Val{2}, dw[2], C, dx))
        return diff + adve
    end
end

function diffusion(box::Box{BT}, d, C, w) where {BT}
    dC = Array{typeof(1.0/u"yr")}(undef, box.grid_size...)
    diffusion!(box, d, C, w, dC)
    return dC
end

function step!(input, state)
    input.solver(
        (a, t) -> diffusion(input.box, input.diffusion_coefficient, a, .-input.topography),
        state.value, state.time, input.time.Δt)
    state.time += input.time.Δt
end

end

module FullTransport

using CarboKitten  #: Box, TimeProperties
using CarboKitten.Stencil: stencil!, Size
using Unitful
using ..DifferentialOperators: central_difference, upwind, laplacian

@kwdef struct Input
    box::Box
    time::TimeProperties
    initial_state::Array{Float64}
    topography::Array{typeof(1.0u"m")}

    diffusion_coefficient
    wave_velocity
    wave_shear

    solver
end

@kwdef mutable struct State
    time::typeof(1.0u"Myr")
    # delta::Array{Float64}
    value::Array{Float64}
end

initial_state(input) = State(
    time = input.time.t0,
    # delta = zeros(Float64, input.box.grid_size...),
    value = copy(input.initial_state))

function transport!(box::Box{BT}, diffusion_coefficient, wave_velocity, wave_shear, C, w, dC) where {BT}
    dx = box.phys_scale
    stencil!(BT, Size(3, 3), dC, C, w) do C, w
        d = diffusion_coefficient
        v = wave_velocity(w[2, 2])
        s = wave_shear(w[2, 2])

        dw = (central_difference(Val{1}, w, dx), central_difference(Val{2}, w, dx))
        diff = -d * C[2, 2] * laplacian(w, dx)
        adve = -d * (upwind(Val{1}, dw[1], C, dx) + upwind(Val{2}, dw[2], C, dx))
        wave_adv = -(upwind(Val{1}, v[1], a, dx) + upwind(Val{2}, v[2], a, dx))
        wave_shr = C[2, 2] * (s[1] * dw[1] + s[2] * dw[2])
        return diff + adve + wave_adve + wave_shr
    end
end

function transport(box::Box{BT}, diffusion_coefficient, wave_velocity, wave_shear, C, w) where {BT}
    dC = Array{typeof(1.0/u"yr")}(undef, box.grid_size...)
    transport!(box, diffusion_coefficient, wave_velocity, wave_shear, C, w, dC)
    return dC
end

function step!(input, state)
    input.solver(
        (a, t) -> transport(
            input.box, input.diffusion_coefficient, input.wave_velocity, input.wave_shear,
            a, .-input.topography),
        state.value, state.time, input.time.Δt)
    state.time += input.time.Δt
end

end

module FlyingCat

using CarboKitten
using FileIO
using GLMakie

using ..Solvers: forward_euler, runge_kutta_4
using ..Advection: Advection
using ..Runner: run_model

const BOX = CarboKitten.Box{Periodic{2}}(grid_size=(256, 256), phys_scale=0.05u"km")
const INPUT = Advection.Input(
    box = BOX,
    time = TimeProperties(Δt=10u"yr", steps=20),
    initial_state = load("data/cat256.pgm")'[:, end:-1:1] .|> Float64,
    wave_velocity = (4u"m/yr", -2u"m/yr"),
    solver = runge_kutta_4(BOX)
)

function run()
    x, y = box_axes(INPUT.box)

    fig = Figure()
    ax1 = Axis(fig[1, 1])
    heatmap!(ax1, x, y, INPUT.initial_state)

    out = run_model(Model{Advection}, INPUT)
    ax2 = Axis(fig[1, 2])
    heatmap!(ax2, x, y, out.value)

    return fig
end

end

module CatTopography

using CarboKitten
using FileIO
using GLMakie

using ..Solvers: forward_euler, runge_kutta_4
using ..Diffusion: Diffusion
using ..Runner: run_model

const BOX = CarboKitten.Box{Periodic{2}}(grid_size=(256, 256), phys_scale=0.05u"km")
const INPUT = Diffusion.Input(
    box = BOX, 
    time = TimeProperties(Δt=10u"yr", steps=20),
    initial_state = ones(Float64, 256, 256),
    topography = (load("data/cat256.pgm")'[:, end:-1:1] .|> Float64) * -50u"m",
    diffusion_coefficient = 20u"m/Myr",
    #wave_velocity = (4u"m/yr", -2u"m/yr"),
    solver = runge_kutta_4(BOX)
)

function run()
    x, y = box_axes(INPUT.box)

    fig = Figure()
    ax1 = Axis(fig[1, 1], aspect=1)
    heatmap!(ax1, x, y, INPUT.topography, colormap=Reverse(:grays))
    out = run_model(Model{Diffusion}, INPUT)
    ax2 = Axis(fig[1, 2], aspect=1)
    heatmap!(ax2, x, y, out.value, colormap=Reverse(:curl))

    return fig
end

end

module ExplodingCat

using CarboKitten
using FileIO
using GLMakie

using ..Solvers: forward_euler, runge_kutta_4
using ..Diffusion: Diffusion
using ..Runner: run_model

const N = 288

function load_cat()
    b = div(N - 256, 2)
    cat = zeros(Float64, N, N)
    cat[b+1:256+b, b+1:256+b] .= (load("data/cat256.pgm")'[:, end:-1:1] .|> Float64)
    return cat
end

const BOX = CarboKitten.Box{Reflected{2}}(grid_size=(N, N), phys_scale=0.05u"km")
const X, Y = box_axes(BOX)
const INPUT = Diffusion.Input(
    box = BOX, 
    time = TimeProperties(Δt=100u"yr", steps=100),
    initial_state = load_cat(),
    topography = ((x, y) -> 30.0u"m" * exp(-((x-6.4u"km")^2 + (y-6.4u"km")^2)/(2*(3.0u"km")^2)) - 30.0u"m").(X, Y'),
    diffusion_coefficient = 30.0u"m/yr",
    #wave_velocity = (4u"m/yr", -2u"m/yr"),
    solver = runge_kutta_4(BOX)
)

function run()
    x, y = box_axes(INPUT.box)

    fig = Figure()
    ax1 = Axis(fig[1, 1], aspect=1)
    hm1 = heatmap!(ax1, x, y, INPUT.topography / u"m")
    Colorbar(fig[2,1], hm1, vertical=false)

    out = run_model(Model{Diffusion}, INPUT)
    ax2 = Axis(fig[1, 2], aspect=1)
    hm2 = heatmap!(ax2, x, y, out.value, colorrange=0.0:2.0)
    Colorbar(fig[2,2], hm2, vertical=false)

    return fig
end

end
# ~/~ end
