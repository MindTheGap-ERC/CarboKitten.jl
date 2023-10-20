# ~/~ begin <<docs/src/ca-with-production.md#examples/cap-extSL.jl>>[init]
using CarboKitten.CaProd
using DataFrames
using CSV
using Interpolations
#datasl = CSV.read("data/bs92-sealevel-curve.csv", DataFrame, header=1)

DEFAULT_INPUT = CaProd.Input(
    sea_level = t -> CSV.read("data/all-sealevel/sl.csv", DataFrame, header=1)[!,1][round(Int,t*1000)+1],
    subsidence_rate=50.0,
    initial_depth=x -> x / 2.0,
    grid_size=(50, 100),
    phys_scale=1,
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

CaProd.main(DEFAULT_INPUT, "data/ca-extSL.h5")
# ~/~ end