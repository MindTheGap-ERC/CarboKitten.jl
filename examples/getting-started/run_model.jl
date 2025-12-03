# ~/~ begin <<docs/src/getting-started.md#examples/getting-started/run_model.jl>>[init]
module Script

using Unitful
using CarboKitten

function main()
    # Configure a list of facies types and their production curves
    facies = [
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

    input = ALCAP.Input(
        # configure a box geometry
        box = Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m"),

        # choose time integration step
        time = TimeProperties(
            Δt=0.0002u"Myr",
            steps=5000),

        # choose what to write to output
        output = Dict(
            # complete in space, sampled in time
            :topography => OutputSpec(write_interval=100),
            # complete in time, but only for a single pixel slice
            :profile    => OutputSpec(slice=(:, 25))
        ),

        # set an initial ramp topography
        initial_topography=(x, y) -> -x / 300.0,

        # set a relative sea level curve
        sea_level = t -> 10.0u"m" * sin(2π * t / 0.20u"Myr") +
                          2.0u"m" * sin(2π * t / 0.03u"Myr"),

        # set a subsidence rate
        subsidence_rate = 50.0u"m/Myr",

        # set the local mean insolation
        insolation = 400.0u"W/m^2",

        # include the facies definition
        facies = facies,

        # set transport properties:
        #   - the sediment buffer keeps track of past sediment
        depositional_resolution = 0.5u"m",
        sediment_buffer_size = 50,

        #   - the disintegration rate sets how fast existing sediment
        #     is removed from the buffer
        disintegration_rate = 50.0u"m/Myr",

        #   - the cementation time sets how fast entrained material
        #     settles, by specifying a half-life time.
        cementation_time = 100u"yr",
    )

    # The following is here for your convenience:
    #   - ensure the output directory exists
    mkpath("data/output")
    #   - do not track the contents of the output directory
    write("data/output/.gitignore", "*\n")

    # run the model
    run_model(Model{ALCAP}, input, "data/output/first-run.h5")
end

end  # of module Script

# Run the script
Script.main()
# ~/~ end
