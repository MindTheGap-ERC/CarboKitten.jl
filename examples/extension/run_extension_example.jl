using CairoMakie
CairoMakie.activate!()
using CarboKitten
import CarboKitten.Output.H5Writer
using CarboKitten.Output.H5Writer: build_environment_grid,write_environment_block!,facies_percent_from_ids, EnvironmentRule,classify_block!,encode_environments, env_selector
using CarboKitten.WaveField: WaveComponent, WaveModel
using CarboKitten.Visualization
using CarboKitten.Visualization: fence_plot, sediment_profile_generic!, extract_column, map_view, wheeler_column, production_curve, summary_plot_from_state, glamour_view, glamour_view!, stratigraphic_column!, resample_to_regular_grid, stratigraphic_column_layers!
using Unitful: ustrip 
using DataFrames
using CarboKitten.Models: ALCAP
using CarboKitten.Boxes: axes as box_axes, Box
using CarboKitten.Export: age_depth_model, DataColumn, read_slice, read_volume, read_header, group_datasets, header_from_input_exact, dataslice_from_state_exact, datavolume_from_state, datacolumn_from_state, data_export, DataSlice
using CarboKitten.Components.WaterDepth: SubsidenceModifier
using CarboKitten.RunExtension: run_single_to_csv, run_once
using GeometryBasics: Vec2
using CSV 
using Interpolations
using LinearAlgebra
using Setfield
using HDF5
using ImageMorphology: label_components
using StatsBase: countmap, mean

function main()
println("Running CarboKitten extension example...")

# --------------------------------------------------
# Paths
# --------------------------------------------------

repo_root = normpath(joinpath(@__DIR__, "..", ".."))
outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

bathymetry_file = joinpath(repo_root, "data", "extension_example", "Bathymetry.csv")
subsidence_file = joinpath(repo_root, "data", "extension_example", "subsidence_map_70x71.csv")
sealevel_file   = joinpath(repo_root, "data", "extension_example", "SeaLevel_0.02Myr.csv")

# --------------------------------------------------
# Domain and time
# --------------------------------------------------

box = Box{Coast}(grid_size = (70,71), phys_scale = 200.0u"m")

time = TimeProperties(
    t0 = -157.3u"Myr",
    Δt = 0.02u"Myr",
    steps = 595
)

# --------------------------------------------------
# Sea level
# --------------------------------------------------

sealevel = CSV.read(sealevel_file, DataFrame)
sort!(sealevel, [:Age_Ma])

sl_time = sealevel.Age_Ma .* u"Myr"
sl_height = sealevel.LT_reconstructed .* u"m"

haq_sea_level = linear_interpolation(sl_time, sl_height, extrapolation_bc = Flat())
sl0 = haq_sea_level(time.t0)
sea_level_shifted = τ -> haq_sea_level(τ) - sl0

# --------------------------------------------------
# Bathymetry and subsidence
# --------------------------------------------------

bathymetry = Matrix(CSV.read(bathymetry_file, DataFrame; header = false))'.* u"m"

subsidence = Matrix(CSV.read(subsidence_file, DataFrame; header = false))
subsidence = Float64.(subsidence) .* u"m/Myr"

@assert size(bathymetry) == box.grid_size
@assert size(subsidence) == box.grid_size
# --------------------------------------------------
# Facies names and colors
# --------------------------------------------------

facies_names = ["Ooids","Calcareous Algaes","Carbonate Mud","Microbialite","Coral-Sponge"]
    FACIES_COLORS  = [:dodgerblue, :green, :brown, :plum, :darkblue]
    facies_cmap = CairoMakie.cgrad(FACIES_COLORS, categorical=true)

# --------------------------------------------------
# Production curves
# --------------------------------------------------

windows_factory1 = []
    windows_factory2 = []
    windows_factory3 = []
    windows_factory4 = []
    windows_factory5 = []

    depths = [0, 5, 10, 20, 50, 100, 150, 200, 300, 400, 500, 600]

    mult1 = [1.0, 1.0, 0.4, 0.05, 0, 0, 0, 0, 0, 0, 0, 0]
    mult2 = [0.8, 1.0, 0.7, 0.2, 0, 0, 0, 0, 0, 0, 0, 0]
    mult3 = [1.0, 1.0, 1.0, 0.7, 0.3, 0.1, 0, 0, 0, 0, 0, 0]
    mult4 = [0.2, 0.5, 0.9, 1.0, 1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    mult5 = [0.0, 0.1, 0.3, 0.6, 1.0, 1.0, 0.9, 0.7, 0.5, 0.3, 0.1, 0.0]

    make_depth_knots(depths, mults) = [(d*u"m", Float64(m)) for (d, m) in zip(depths, mults)]
    make_time_windows(w) = [(t1*u"Myr", t2*u"Myr", Float64(m)) for (t1, t2, m) in w]

# --------------------------------------------------
# Output specification
# --------------------------------------------------

output_spec = Dict(
    :facies => OutputSpec(write_interval = 200),
    :topk => OutputSpec(write_interval = 200),
    :topography => OutputSpec(write_interval = 100),
    :profile => OutputSpec(slice = (:, 25), write_interval = 1),
)

# --------------------------------------------------
# Wave model
# --------------------------------------------------

H1 = 4.45u"m"
A1 = H1 / 2
az1 = 245

H2 = 2u"m"
A2 = H2 / 2
az2 = 60

comp1 = WaveComponent(
    A1,
    12.12u"s",
    deg2rad(90 - az1),
    0.0,
    0.0u"1/m"
)

comp2 = WaveComponent(
    A2,
    10u"s",
    deg2rad(90 - az2),
    0.0,
    0.0u"1/m"
)

wm = WaveModel([comp1, comp2], 1.0)

# --------------------------------------------------
# Facies
# --------------------------------------------------
	run_id=1
	phi0_1 = 0.45
	phi0_2 = 0.5
	phi0_3 = 0.55
	phi0_4 = 0.7
	phi0_5 = 0.45
	facies = [
        ALCAP.Facies(name="ooids",   
					 maximum_growth_rate=120.0u"m/Myr", 
					 depth_knots=make_depth_knots(depths, mult1), 
					 time_windows=make_time_windows(windows_factory1),
          		     diffusion_coefficient=200.0u"m/Myr", 
					 neighbour_radius=3, 
            		 activation_range=(12,30),
            		 viability_range=(8,40),
					 active = true,
					 depositional_porosity = phi0_1,
					compaction_curve = z -> phi0_1 * exp(-0.4 * z)
					),
        ALCAP.Facies(name="calcareous algeas",  
					 maximum_growth_rate=120.0u"m/Myr",
					 depth_knots=make_depth_knots(depths, mult2), time_windows=make_time_windows(windows_factory2),
					diffusion_coefficient=25.0u"m/Myr", 
					 neighbour_radius=2,
					activation_range=(4,10), 
					viability_range=(3,20),
					 active = true,
					depositional_porosity = phi0_2,
					compaction_curve = z -> phi0_2 * exp(-0.45 * z)
					),
        ALCAP.Facies(name="carbonate mud",    
					 maximum_growth_rate=80.0u"m/Myr",
					 depth_knots=make_depth_knots(depths, mult3), time_windows=make_time_windows(windows_factory3), diffusion_coefficient=300.0u"m/Myr",
					 neighbour_radius=3,
					 viability_range=(2,48), 
					 activation_range=(2,48), 
					 active=true,
					depositional_porosity = phi0_3,
					compaction_curve = z -> phi0_3 * exp(-0.4 * z)
					),
        ALCAP.Facies(name="microbialites",  
					 maximum_growth_rate=80.0u"m/Myr",
					 depth_knots=make_depth_knots(depths, mult4), time_windows=make_time_windows(windows_factory4),
					 diffusion_coefficient=50.0u"m/Myr", 
					 neighbour_radius=3,
					 viability_range=(2,16), 
					 activation_range=(3,12), 
					 active=true,
					depositional_porosity = phi0_4,
					compaction_curve = z -> phi0_4 * exp(-0.7 * z)
					),
        ALCAP.Facies(name="coral-sponge", 
					 maximum_growth_rate=30.0u"m/Myr",
					 depth_knots=make_depth_knots(depths, mult5), time_windows=make_time_windows(windows_factory5),
            		diffusion_coefficient=25.0u"m/Myr",
					 neighbour_radius=2,
					 activation_range=(4,10), 
					viability_range=(3,10), 
					 active=true,
					depositional_porosity = phi0_5,
					compaction_curve = z -> phi0_5 * exp(-0.8 * z)
					),
    ]

# --------------------------------------------------
# Input
# --------------------------------------------------

input = ALCAP.Input(
        tag = "morris_$(run_id)",
        time = time,
        box = box,
        facies = facies,
        ca_interval = 5,
        ca_random_seed = 1,
        sea_level = sea_level_shifted,
        initial_topography = bathymetry,
        subsidence_rate = subsidence,
	subsidence_modifiers = [],   
	 wave_model = wm,
        insolation = 400.0u"W/m^2",
        sediment_buffer_size = 50,
        depositional_resolution = 1.5u"m",
        disintegration_rate =10u"m/Myr",
        cementation_time = 0.2u"Myr",
        intertidal_zone = 1u"m",
        output = output_spec
    )

# --------------------------------------------------
# Run model
# --------------------------------------------------

env_names = String[
    "Ooid Shoal",
    "Algal Reef",
    "Microbialite",
    "Coral-Sponge Mound",
    "Open Lagoon",
    "Restricted Lagoon",
    "Semi-Closed Lagoon",
    "Unclassified"
]
	out_csv  = joinpath(outdir, "example_run_metrics.csv")
	fail_csv = joinpath(outdir, "example_run_failures.csv")



n_facies = length(facies)

out = run_single_to_csv(
    input;
    facies_names = facies_names,
    env_names = env_names,
    out_csv = out_csv,
    fail_csv = fail_csv,
    run_id = run_id,
    n_facies = n_facies
)
	state = out.state

# --------------------------------------------------
# Facies plots
# --------------------------------------------------

Nx, Ny, _ = size(state.block_cube)
xs = round.(Int, range(1, Nx; length=6))
ys = round.(Int, range(1, Ny; length=6))

grid = resample_to_regular_grid(
    state.layer_facies_hist,
    state.layer_thickness_hist,
    1.0f0
)

fig_fence_facies = fence_plot(
    grid;
    category_names = facies_names,
	facies_colors = FACIES_COLORS,
    dx = 200.0,
    dz = 1.0
)

save(joinpath(outdir, "fence_facies.png"), fig_fence_facies)

zm = 250  # meters
dz_m = ustrip(input.depositional_resolution / u"m")	
fig_map_facies = map_view(
    grid;
    z0 = Int(round(zm / dz_m)),
    category_names = facies_names,
    category_colors = FACIES_COLORS,
    colorbar_label = "Facies"
)
save(joinpath(outdir, "map_facies.png"), fig_map_facies)

# --------------------------------------------------
# Environment classification
# --------------------------------------------------

env_rules = [
    EnvironmentRule(
        "Ooid Shoal",
        Dict(1 => (70.0, 100.0)),
        nothing,
        (180.0717, 215.471)
    ),
    EnvironmentRule(
        "Algal Reef",
        Dict(2 => (40.0, 100.0)),
        nothing,
        nothing
    ),
    EnvironmentRule(
        "Microbialite",
        Dict(4 => (60.0, 100.0)),
        nothing,
        nothing
    ),
    EnvironmentRule(
        "Coral-Sponge Mound",
        Dict(5 => (50.0, 100.0)),
        nothing,
        nothing
    ),
    EnvironmentRule(
        "Open-Lagoon",
        Dict(),
        (24, 127.4029),
        nothing
    ),
    EnvironmentRule(
        "Restricted Lagoon",
        Dict(),
        (0, 8),
        nothing
    ),
    EnvironmentRule(
        "Semi-Closed Lagoon",
        Dict(),
        (8, 24),
        nothing
    ),
]

ENV_COLORS  = [:dodgerblue, :coral, :darkkhaki, :purple,:yellow, :green, :goldenrod, :lightgray]
env_grid = build_environment_grid(
    state.layer_facies_hist,
    state.layer_thickness_hist,
    state.layer_wdepth_hist,
    state.layer_energy_hist;
    Δi = 5,
    Δj = 5,
    Δk = 5,
    env_rules = env_rules,
    n_facies = length(input.facies)
)


encoded = encode_environments(env_grid, env_names)

# --------------------------------------------------
# Environment plots
# --------------------------------------------------

fig_fence_env = fence_plot(
    encoded;
    category_names = env_names,
	facies_colors = ENV_COLORS,
    dx = 200.0,
    dz = 1.5,
	VE = 1
)
save(joinpath(outdir, "fence_environment.png"), fig_fence_env)

fig_map_env = map_view(
    encoded;
    z0 = Int(round(zm / dz_m)),
    category_names = env_names,
    category_colors = ENV_COLORS,
    colorbar_label = "Environments"
)
save(joinpath(outdir, "map_environment.png"), fig_map_env)

println("Done.")
println("Outputs written to: ", outdir)
end

main()