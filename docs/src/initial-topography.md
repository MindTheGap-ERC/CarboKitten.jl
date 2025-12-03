# Generate Initial topography

Many carbonate platforms grow on a pre-exisiting abandoned carbonate platform. Therefore, it is reasonable to start the simulation with a platform-shaped initial topography.

The workflow to generate the topography is:

1. run the CarboKitten.jl to get the results.
2. export the sedimentation for each cell.
3. calculate the elvelation for each cell by adding subsidence. 

``` {.julia file=examples/initial_topography/get_init_topo.jl}
module Prerun

using CarboKitten
using CarboKitten.Export: read_volume

using Tables
using CSV: write as write_csv

const PATH = "data/output/initial-topography"
const TAG = "prerun"
const DATAFILE = joinpath(PATH, "initial-topography.csv")

<<example-init-topo>>
end  # module Prerun
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
    tag=TAG,
    box=Box{Coast}(grid_size=(100, 70), phys_scale=170u"m"),
    time=TimeProperties(
        Δt=0.0001u"Myr",
        steps=5000
    ),
    output=Dict(
        :topography => OutputSpec(write_interval=1000)
    ),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level= t -> 0.0u"m",
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    facies=FACIES)

```

For the second step, the disintegration, production, deposition and sedimentation are exported respectively. The starting bathymetry for this run is set to be a slight slope (slope: 1/300)

``` {.julia #example-init-topo}
function prerun()
    mkpath(PATH)
    run_model(Model{ALCAP}, INPUT, "$(PATH)/$(TAG).h5")
end

function save_final_topography(prerun_filename)
    header, data = read_volume(prerun_filename, :topography)

    t = header.axes.t
    h0 = header.initial_topography
    subsidence = header.subsidence_rate * (t[end] - t[1])
    delta_h = data.sediment_thickness[:, :, end]
    h = h0 .+ delta_h .- subsidence

    write_csv(DATAFILE, h |> in_units_of(u"m") |> Tables.table)
end
```

For the third step, the elevation is calculated through by substracting the sedimentation with subsidence. 

``` {.julia #example-init-topo}
function calculate_bathymetry(data,INPUT)
    Bathy = zeros(x_grid_size, y_grid_size) .*1.0u"m"
    Bathy .= starting_bathy() .+ data .- INPUT.subsidence_rate .* INPUT.time.Δt .* INPUT.time.steps
    OfficialCSV.write("examples/initial_topography/example_init_topo.csv", DataFrame(Bathy,:auto))
end

calculate_bathymetry(data_sed,INPUT)
```

The resultant initial topography that's ready for your run is stored in csv format. 

The next step is to import the csv file, through `initial_topography(path)`, where you may need to secify the location of the csv file. 

