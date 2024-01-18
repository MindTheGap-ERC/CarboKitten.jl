# ~/~ begin <<docs/src/ca-with-production.md#examples/caps-osc.jl>>[init]
using CarboKitten.CaProd
using CarboKitten.Burgess2013
using CSV
using DataFrames
using Interpolations

function sealevel_curve(t,filepath)
    data = DataFrame(CSV.File(filepath))
    data = hcat(collect(0:length(data[:,1]).-1)./1000, data[:,1]) #The output sealevel curve from R does not have time tab and have to add it 
    x = linear_interpolation(data[:,1], data[:,2])
    return x(t)
end 


DEFAULT_INPUT = CaProd.Input(
    sea_level = t -> sealevel_curve(t,"data/miller.csv"),# 
    subsidence_rate = 50.0,
    initial_depth = x -> x / 2,
    grid_size = (50, 100),
    phys_scale = 1.0,
    Δt = 0.001,
    write_interval = 1,
    time_steps = 1000,
    facies = [
        Burgess2013.Facies((4, 10), (6, 10), 500.0, 0.8, 300, 1000, 2730, 0.5),
        Burgess2013.Facies((4, 10), (6, 10), 400.0, 0.1, 300, 1000, 2730, 0.5),
        Burgess2013.Facies((4, 10), (6, 10), 100.0, 0.005, 300, 1000, 2730, 0.5)
    ],
    insolation = 2000.0,
    temp = 288.0,
    precip = 1000.0,
    pco2 = 10^(-1.5),
    alpha = 2e-3

)

CaProd.main(DEFAULT_INPUT, "data/caps-test.h5")
# ~/~ end