# ~/~ begin <<docs/src/ca-with-production.md#examples/production-only/cap-slope.jl>>[init]
#| creates: data/ca-prod-slope.h5
#| requires: src/CaProd.jl

using CarboKitten.CaProd

DEFAULT_INPUT = CaProd.Input(
    sea_level=_ -> 0.0,
    subsidence_rate=50.0,
    initial_depth=x -> x / 2.0,
    grid_size=(50, 100),
    phys_scale=1.0,
    Δt=0.001,
    write_interval=1,
    time_steps=1000,
    facies=[
        CaProd.Facies((4, 10), (6, 10), 500.0, 0.8, 300),
        CaProd.Facies((4, 10), (6, 10), 400.0, 0.1, 300),
        CaProd.Facies((4, 10), (6, 10), 100.0, 0.005, 300)
    ],
    insolation=2000.0
)

CaProd.main(DEFAULT_INPUT, "data/ca-prod-slope.h5")
# ~/~ end