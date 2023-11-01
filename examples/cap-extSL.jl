# ~/~ begin <<docs/src/ca-with-production.md#examples/cap-extSL.jl>>[init]
using CarboKitten.CaProd
using DataFrames
using CSV
using Interpolations


function sealevel_curve(t,filepath)
    data = DataFrame(CSV.File(filepath))
    data = hcat(collect(0:length(data.sl).-1)./1000, data.sl) #The output sealevel curve from R does not have time tab and have to add it 
    x = linear_interpolation(data[:,1], data[:,2])
    return x(t)
end 

DEFAULT_INPUT = CaProd.Input(
    sea_level = t -> sealevel_curve(t,file_path),
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

folder_path = "data/all-sealevel"
files = filter(x -> occursin(".csv", x), readdir(folder_path))

for i in files
    file_path = joinpath(folder_path,i)
    name, _ = splitext(i)
    CaProd.main(DEFAULT_INPUT, "data/ca-extSL-$name.h5")
end
