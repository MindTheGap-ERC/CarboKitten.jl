# ~/~ begin <<docs/src/finite-difference-transport.md#examples/transport/fd.jl>>[init]
using Unitful
using CarboKitten



module CatTopography

using CarboKitten
using FileIO
using GLMakie

using ..Solvers: forward_euler, runge_kutta_4
using ..Diffusion: Diffusion
using ..Transport: Transport
using ..Runner: run_model

const BOX = CarboKitten.Box{Periodic{2}}(grid_size=(256, 256), phys_scale=0.05u"km")
const INPUT = Transport.Input(
    box = BOX, 
    time = TimeProperties(Δt=10u"yr", steps=20),
    initial_state = ones(Float64, 256, 256),
    topography = (load("data/cat256.pgm")'[:, end:-1:1] .|> Float64) * -50u"m",
    diffusivity = 20u"m/Myr",
    #wave_velocity = (4u"m/yr", -2u"m/yr"),
    solver = runge_kutta_4(BOX)
)

function run()
    x, y = box_axes(INPUT.box)

    fig = Figure()
    ax1 = Axis(fig[1, 1], aspect=1)
    heatmap!(ax1, x, y, INPUT.topography, colormap=Reverse(:grays))
    out = run_model(Model{Transport}, INPUT)
    ax2 = Axis(fig[1, 2], aspect=1)
    heatmap!(ax2, x, y, out.value, colormap=Reverse(:curl))

    return fig
end

end
# ~/~ end
