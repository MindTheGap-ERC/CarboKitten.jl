# ~/~ begin <<docs/src/initial-topography.md#examples/initial_topography/get_init_topo.jl>>[init]
module Prerun

using CarboKitten
using CarboKitten.Export: read_volume

using Tables
using CSV: write as write_csv

const PATH = "data/output/initial-topography"
const TAG = "prerun"
const DATAFILE = joinpath(PATH, "initial-topography.csv")

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

# ~/~ end
# ~/~ begin <<docs/src/initial-topography.md#example-init-topo>>[1]
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
# ~/~ end
end  # module Prerun
# ~/~ end
