# ~/~ begin <<docs/src/initial-topography.md#examples/initial_topography/production_run.jl>>[init]
module ProductionRun

using CarboKitten
using CSV
using Tables

include("get_init_topo.jl")

const PATH = Prerun.PATH

# ~/~ begin <<docs/src/initial-topography.md#load-initial-topography>>[init]
function load_initial_topography()
    (CSV.File(Prerun.DATAFILE) |> Tables.matrix) * u"m"
end
# ~/~ end

function run()
    input = ALCAP.Input(
        tag = "mainrun",
        time = time = TimeProperties(
            Δt = 200u"yr",
            steps = 5000
        ),
        box = Prerun.INPUT.box,
        facies = Prerun.INPUT.facies,
        sea_level = t -> 0.0u"m",
        initial_topography = load_initial_topography(),
        output=Dict(
            :topography => OutputSpec(write_interval=100),
            :profile => OutputSpec(slice = (:, 35))
        ),
        subsidence_rate = 50.0u"m/Myr",
        insolation = 400.0u"W/m^2",

        transport_solver = Val{:forward_euler},
        sediment_buffer_size = 50,
        depositional_resolution = 0.5u"m",
        cementation_time = 50.0u"yr",
        disintegration_rate = 100.0u"m/Myr")

    run_model(Model{ALCAP}, input, joinpath(PATH, "mainrun.h5"))
end

end  # module MainRun
# ~/~ end
