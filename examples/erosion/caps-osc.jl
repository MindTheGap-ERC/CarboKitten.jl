# ~/~ begin <<docs/src/ca-prod-with-erosion.md#examples/erosion/caps-osc.jl>>[init]
using CarboKitten
using CarboKitten.CaProdErosion
using CarboKitten.Burgess2013

#change sinouid function
DEFAULT_INPUT = CarboKitten.CaProdErosion.Input(
    sea_level = t -> 20 * sin(2π * t) + 5 * sin(2π * t / 0.112), 
    subsidence_rate = 50.0,
    initial_depth = x -> x / 2,
    grid_size = (50, 100),
    phys_scale = 1.0,
    Δt = 0.001,
    write_interval = 1,
    time_steps = 2000,
    facies = [
        Burgess2013.Facies((4, 10), (6, 10), 500.0, 0.8, 300, 1000, 2730, 0.5),
        Burgess2013.Facies((4, 10), (6, 10), 400.0, 0.1, 300, 1000, 2730, 0.5),
        Burgess2013.Facies((4, 10), (6, 10), 100.0, 0.005, 300, 1000, 2730, 0.5)
    ],
    insolation = 2000.0,
    temp = 288.0,
    precip = 1000.0,
    pco2 = 10^(-1.5),
    alpha = 2e-3,
    erosion_type = 1

)

CarboKitten.CaProdErosion.main(DEFAULT_INPUT, "data/capes-osc.h5")
# ~/~ end