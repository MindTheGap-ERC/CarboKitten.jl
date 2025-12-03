# ~/~ begin <<docs/src/initial-topography.md#examples/initial_topography/get_init_topo.jl>>[init]
using CarboKitten

using Unitful
using CarboKitten.Export: data_export, CSV as CSVCarbo
using CarboKitten.Export: read_volume
using HDF5
using DataFrames
import CSV as OfficialCSV

const PATH = "examples/initial_topography"
const TAG = "example_init_topo"
const x_grid_size = 100
const y_grid_size = 70
const PHYS_SCALE = 170.0u"m"
const SLOPE = 300.0
# ~/~ begin <<docs/src/initial-topography.md#example-init-topo>>[init]
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
    box=Box{Coast}(grid_size=(x_grid_size, y_grid_size), phys_scale=PHYS_SCALE),
    time=TimeProperties(
        Δt=0.0001u"Myr",
        steps=5000),
    ca_interval=1,
    initial_topography=(x, y) -> -x / SLOPE,
    sea_level= t -> 0.0u"m",
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    facies=FACIES)

run_model(Model{ALCAP}, INPUT, "examples/initial_topography/example_init_topo.h5")

# ~/~ end
# ~/~ begin <<docs/src/initial-topography.md#example-init-topo>>[1]
function extract_topography()
        data = read_volume("examples/initial_topography/example_init_topo.h5", :full)[2]

        sediment_height = data.sediment_thickness[:, :, end]

        data_sed = DataFrame(
           sediment_height, :auto
        )
        return data_sed
end


data_sed = extract_topography()

function starting_bathy()
    init = ones(x_grid_size, y_grid_size) .*1.0u"m"
    for i in CartesianIndices(init)
        init[i]   = -(i[1]-1) .* PHYS_SCALE ./ SLOPE
    end
    return init
end
# ~/~ end
# ~/~ begin <<docs/src/initial-topography.md#example-init-topo>>[2]
function calculate_bathymetry(data,INPUT)
    Bathy = zeros(x_grid_size, y_grid_size) .*1.0u"m"
    Bathy .= starting_bathy() .+ data .- INPUT.subsidence_rate .* INPUT.time.Δt .* INPUT.time.steps
    OfficialCSV.write("examples/initial_topography/example_init_topo.csv", DataFrame(Bathy,:auto))
end

calculate_bathymetry(data_sed,INPUT)
# ~/~ end
# ~/~ end
