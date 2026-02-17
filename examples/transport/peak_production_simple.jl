module PeakProductionSimple

include("CustomProductionModel.jl")

using CarboKitten
using CarboKitten: Box, TimeProperties, OutputSpec, Model, run_model
using CarboKitten.BoundaryTrait: Reflected
using Unitful
using .CustomProduction: CustomProduction as M

function run_peak_diffusion(;
        dt,
        diffusivity,
        disintegration_rate,
        cementation_time,
        patch_width=2.0u"km",
        production_steps=10,
        t_end=0.1u"Myr",
        output_file="peak_diffusion.h5"
    )

    facies = [M.Facies(diffusion_coefficient=diffusivity)]

    box = Box{Reflected{2}}(
        grid_size=(500, 1), phys_scale=30.0u"m")
    time = TimeProperties(
        Δt=dt,
        steps=(t_end / dt) |> round |> Int)

    centre = box.grid_size[1] * box.phys_scale / 2.0

    # Production: rectangular patch, active only for first `production_steps` steps
    production(step, x, y, _w) =
        (step < production_steps && abs(x - centre) < patch_width) ?
            100.0u"m/Myr" * time.Δt :
            0.0u"m"

    write_interval = div(time.steps, 100)

    input = M.Input(
        box=box,
        time=time,
        output=Dict(
            :profile => OutputSpec(slice=(:, 1), write_interval=write_interval)),
        initial_topography=(_, _) -> -100.0u"m",
        sea_level=t -> 0.0u"m",
        subsidence_rate=0.0u"m/Myr",
        disintegration_rate=disintegration_rate,
        sediment_buffer_size=50,
        depositional_resolution=0.5u"m",
        cementation_time=cementation_time,
        transport_solver=Val{:forward_euler},
        save_active_layer=false,
        facies=facies,
        production=production)

    run_model(Model{M}, input, output_file)
end

function example_1d(;
        dt=0.001u"Myr",
        diffusivity=10.0u"m/yr",
        disintegration_rate=0.01u"m/Myr",
        cementation_time=0.05u"Myr",
        production_steps=10,
        t_end=0.5u"Myr",
        output_file="data/output/example_1d_peak.h5"
    )
    run_peak_diffusion(;
        dt, diffusivity, disintegration_rate, cementation_time,
        production_steps, t_end, output_file
    )
end

end # module
