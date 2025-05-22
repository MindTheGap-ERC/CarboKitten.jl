using Test
using CarboKitten
using Unitful
using CarboKitten.Components
using CarboKitten.Components.Common
using CarboKitten.Components.Denudation
using CarboKitten.Models: WithDenudation as WDn
using CarboKitten.Export: data_export, CSV as CSVCarbo
using CarboKitten.Denudation
using DataFrames
import CSV as Official_CSV
using HDF5: h5open
using GLMakie
const m = u"m"
const Myr = u"Myr"

const PATH = "data/output"
const TAG = "dissolution"

const FACIES = [
    WDn.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=10000u"m",
        reactive_surface=10u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.23u"m/yr"
        ),
]

const PERIOD = 0.2Myr
const AMPLITUDE = 20.0m
const DENUDATION = Dissolution(temp=293.0u"K", precip=1.0u"m", pco2=10^(-2.5) * u"atm", reactionrate=2e-3u"m/yr")


const INPUT = WDn.Input(
    tag="$TAG",
    box=Box{Coast}(grid_size=(100, 50), phys_scale=150.0m),
    time=TimeProperties(
        Δt=0.0002Myr,
        steps=1000,
        write_interval=1),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=t -> AMPLITUDE * sin(2π * t / PERIOD),
    subsidence_rate=20.0m / Myr,
    disintegration_rate=500.0m / Myr,
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5m,
    facies=FACIES,
    denudation = DENUDATION)

function main()
    H5Writer.run_model(Model{WDn}, INPUT, "$(PATH)/$(TAG).h5")

    data_export(
        CSVCarbo(tuple.(10:10:100, 25),
          :sediment_accumulation_curve => "$(PATH)/$(TAG)_sac.csv",
          :age_depth_model => "$(PATH)/$(TAG)_adm.csv",
          :stratigraphic_column => "$(PATH)/$(TAG)_sc.csv",
          :metadata => "$(PATH)/$(TAG).toml"),
        "$(PATH)/$(TAG).h5")
end

main()

function read_csvfiles(PATH::String,TAG::String)
    sc = Official_CSV.read("$(PATH)/$(TAG)_sc.csv", DataFrame)
    sac = Official_CSV.read("$(PATH)/$(TAG)_sac.csv", DataFrame)
    return sc, sac
end

sc,sac = read_csvfiles(PATH,TAG)

@test sum(sc.var"sc1_f1 [m]") ≈ sac.var"sac1 [m]"[end]


## find the largest discrepancy in the y-axis
function calculate_discrepancy(sc_col, sac_col)
    sum(sc_col) ./sac_col[end]
end

function calculate_discrepancy_all(sc::DataFrame,sac::DataFrame)
    discrepancy_matrix = Float64[]
    for i in 2:size(sc)[2]
        a = calculate_discrepancy(sc[:,i], sac[:,i])
        push!(discrepancy_matrix, a)
    end
    return discrepancy_matrix
end

function find_discrepancy_location(SC::DataFrame,SAC::DataFrame)
    discrepancy_matrix = calculate_discrepancy_all(sc, sac)
    @show discrepancy_matrix
    discrep_percent,location = findmax(discrepancy_matrix)
    return (discrep_percent, location)
end

find_discrepancy_location(sc,sac)

## extract deposition, production and disintegration with time from HDF5

function extract_data_from_hdf5(FILE::String)
    h5open(FILE, "r") do fid
        production = read(fid["/production"], Float64)[1,100,25,:]
        deposition = read(fid["/deposition"], Float64)[1,100,25,:]
        disintegration = read(fid["/disintegration"], Float64)[1,100,25,:]
    
        @show size(production)
    data_prod = DataFrame(
        production = production,
    )

    data_depo = DataFrame(
        deposition = deposition,
    )   

    data_dis = DataFrame(
        disintegration = disintegration,
    )

    return data_prod, data_depo, data_dis
    end

end

data_prod,data_depo,data_dis = extract_data_from_hdf5("$(PATH)/$(TAG).h5") 

@test sum(data_depo.deposition) ≈ sum(data_prod.production) - sum(data_dis.disintegration)

## Plotting production, deposition and disintegration
fig = Figure()
ax0 = Axis(fig[1, 1], xlabel="Time step", ylabel="Production [m]", title="Production")
ax2 = Axis(fig[2, 1], xlabel="Time step", ylabel="Deposition [m]", title="Deposition")
ax3 = Axis(fig[3, 1], xlabel="Time step", ylabel="Disintegration [m]", title="Disintegration")
lines!(ax0,data_prod.production, label="Production", color=:blue)
lines!(ax2,data_depo.deposition, label="Deposition", color=:red) 
lines!(ax3,data_dis.disintegration, label="Disintegration", color=:green)
save("$(PATH)/$(TAG).png",fig)
