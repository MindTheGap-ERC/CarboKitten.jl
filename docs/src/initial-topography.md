# Generate Initial topography

Many carbonate platforms grow on a pre-exisiting abandoned carbonate platform. Therefore, it is reasonable to start the simulation with a platform-shaped initial topography.

The workflow to generate the topography is:

1. run the CarboKitten.jl to get the results.
2. export the sedimentation for each cell.
3. calculate the elvelation for each cell by adding subsidence. 

``` {.julia file=examples/initial_topography/get_init_topo.jl}
using CarboKitten

using Unitful
using CarboKitten.Export: data_export, CSV as CSVCarbo
using HDF5
using DataFrames
import CSV as OfficialCSV

const PATH = "data/init_topo"

const TAG = "example_init_topo"


<<example-init-topo>>
```

In the first step, the code is listed below. In this example run, constant sea-level has been used. One thing should be noticed is that the grids are 100 by 70, with scale of 170 m. These values should be same as your runs later.

``` {.julia #example-init-topo}
const FACIES = [
    ALCAP.Facies(
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=50.0u"m/yr"),
    ALCAP.Facies(
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=25.0u"m/yr"),
    ALCAP.Facies(
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
```

For the second step, the disintegration, production, deposition and sedimentation are exported respectively. The starting bathymetry for this run is set to be a slight slope (slope: 1/300)

``` {.julia #example-init-topo}
function extract_topography(path,tag)
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
```

For the third step, the elevation is calculated through by substracting the sedimentation with subsidence. 

``` {.julia #example-init-topo}
function calculate_bathymetry(data,INPUT)
    Bathy = zeros(100, 70) .*1.0u"m"
    Bathy .= starting_bathy() .+ data .- INPUT.subsidence_rate .* INPUT.time.Δt .* INPUT.time.steps
    OfficialCSV.write("$(PATH)/$(TAG).csv", DataFrame(Bathy,:auto))
end
```

The resultant initial topography that's ready for your run is stored in csv format. 

The next step is to import the csv file, through `initial_topography(path)`, where you may need to secify the location of the csv file. 

