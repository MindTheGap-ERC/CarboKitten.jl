# ~/~ begin <<docs/src/finite-difference-transport.md#examples/transport/flying_cat.jl>>[init]
module FlyingCat

include("runner.jl")
include("test_model.jl")

using CarboKitten
using FileIO
using GLMakie

using CarboKitten.Transport.Solvers: forward_euler, runge_kutta_4

const BOX = CarboKitten.Box{Periodic{2}}(
    grid_size=(256, 256), phys_scale=0.05u"km")

const INPUT = TestModel.Input(
    box = BOX,
    time = TimeProperties(Δt=10u"yr", steps=50),
    topography = zeros(typeof(1.0u"m"), BOX.grid_size),
    initial_state = load("data/cat256.pgm")'[:, end:-1:1] .|> Float64,
    wave_velocity = _ -> (4u"m/yr", -3.0u"m/yr"),
    solver = runge_kutta_4(BOX)
)

function run()
    x, y = box_axes(INPUT.box)

    fig = Figure()
    ax1 = Axis(fig[1, 1], aspect=1)
    hm1 = heatmap!(ax1, x, y, INPUT.initial_state, colorrange=(0.0,0.7))
    Colorbar(fig[2, 1], hm1, vertical=false)

    out = Runner.run_model(Model{TestModel}, INPUT)
    ax2 = Axis(fig[1, 2], aspect=1)
    hm2 = heatmap!(ax2, x, y, out.value, colorrange=(0.0,0.7))
    Colorbar(fig[2, 2], hm2, vertical=false)

    return fig
end

end
# ~/~ end
