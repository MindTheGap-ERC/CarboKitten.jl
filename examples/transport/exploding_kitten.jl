# ~/~ begin <<docs/src/algorithms/finite-difference-transport.md#examples/transport/exploding_kitten.jl>>[init]

module ExplodingKitten

include("test_model.jl")

using CarboKitten
using FileIO
using GLMakie

using CarboKitten.Transport.Solvers: runge_kutta_4, forward_euler

GLMakie.activate!()

const N = 288

function load_cat()
    b = div(N - 256, 2)
    cat = zeros(Float64, N, N)
    cat[b+1:256+b, b+1:256+b] .= (load("data/cat256.pgm")'[:, end:-1:1] .|> Float64)
    return cat
end

const BOX = CarboKitten.Box{Reflected{2}}(grid_size=(N, N), phys_scale=0.05u"km")
const X, Y = box_axes(BOX)
const INPUT = TestModel.Input(
    box = BOX, 
    time = TimeProperties(Δt=100u"yr", steps=100),
    initial_state = load_cat(),
    topography = ((x, y) -> 30.0u"m" * exp(-((x-7.2u"km")^2 + (y-7.2u"km")^2)/(2*(3.0u"km")^2)) - 30.0u"m").(X, Y'),
    diffusivity = 30.0u"m/yr",
    solver = forward_euler
)

function run()
    x, y = box_axes(INPUT.box)

    fig = Figure(size=(800, 400))
    ax1 = Axis(fig[1, 1], aspect=1)
    hm1 = heatmap!(ax1,
        x |> in_units_of(u"km"), y |> in_units_of(u"km"),
        INPUT.topography / u"m")
    Colorbar(fig[2,1], hm1, vertical=false)

    do_nothing(_, _) = nothing
    out = run_model(do_nothing, Model{TestModel}, INPUT)
    ax2 = Axis(fig[1, 2], aspect=1)
    hm2 = heatmap!(ax2,
        x |> in_units_of(u"km"), y |> in_units_of(u"km"),
        out.value, colorrange=(0.0,0.8))
    Colorbar(fig[2,2], hm2, vertical=false)

    save("docs/src/_fig/exploding_kitten.png", fig)
end

end

ExplodingKitten.run()
# ~/~ end
