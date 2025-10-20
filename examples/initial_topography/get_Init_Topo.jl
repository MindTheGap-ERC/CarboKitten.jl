using CarboKitten

using Unitful
using CarboKitten.Export: data_export, CSV as CSVCarbo
using HDF5
using DataFrames
import CSV as OfficialCSV

const PATH = "data/init_topo"

const TAG = "example_init_topo"

const FACIES = [
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=50.0u"m/yr"),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=25.0u"m/yr"),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=12.5u"m/yr")
]

const INPUT = ALCAP.Input(
    tag="$TAG",
    box=Box{Coast}(grid_size=(100, 70), phys_scale=170.0u"m"),
    time=TimeProperties(
        Δt=0.0001u"Myr",
        steps=5000),
    #output=Dict(
    #        :profile => OutputSpec(slice = (:, 35), write_interval = 1)),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level= t -> 0.0u"m",
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    facies=FACIES)

run_model(Model{ALCAP}, INPUT, "$(PATH)/$(TAG).h5")

function extract_topography(PATH,TAG)
    h5open("$(PATH)/$(TAG).h5", "r") do fid
        disintegration = read(fid["full/disintegration"])[1,:,:,end]
        @show size(disintegration)
        production = read(fid["full/production"])[1,:,:,end]
        deposition = read(fid["full/deposition"])[1,:,:,end]
        sediment_height = read(fid["full/sediment_thickness"])[:,:,end]
        @show size(sediment_height)
        data_dis = DataFrame(
            disintegration, :auto
        )
        @show size(data_dis)
        data_pro = DataFrame(
            production, :auto
        )   
        data_dep = DataFrame(
           deposition, :auto
        )
        data_sed = DataFrame(
           sediment_height, :auto
        )
        return data_dis.*1.0u"m", data_pro.*1.0u"m", data_dep.*1.0u"m", data_sed.*1.0u"m"
end
end

data_dis, data_pro, data_dep, data_sed = extract_topography(PATH,TAG)

function starting_bathy()
    init = ones(100, 70) .*1.0u"m"
    for i in CartesianIndices(init)
        init[i]   = -(i[1]-1) .* 170u"m" ./ 300
    end
    return init
end

starting_bathy()

function calculate_bathymetry(data,INPUT)
    Bathy = zeros(100, 70) .*1.0u"m"
    Bathy .= starting_bathy() .+ data .- INPUT.subsidence_rate .* INPUT.time.Δt .* INPUT.time.steps
    OfficialCSV.write("$(PATH)/$(TAG).csv", DataFrame(Bathy,:auto))
end


calculate_bathymetry(data_sed,INPUT)



