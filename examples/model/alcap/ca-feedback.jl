# ~/~ begin <<docs/src/components/cellular-automata.md#examples/model/alcap/ca-feedback.jl>>[init]
module Script
    using CarboKitten
    using CarboKitten.Production
    using CarboKitten.Models: ALCAP as M

    initial_topography(x, y) =
        - sqrt((x - 7.5u"km")^2 + (y - 7.5u"km")^2) / 100.0

    function main()
        res = 200
        steps = 5000
        phys_scale = 15.0u"km" / res

        output = Dict(
            :topography => OutputSpec(write_interval = max(1, div(steps, 50))),
            :profile    => OutputSpec(slice = (:, div(res, 2)+1)))

        facies = [
            M.Facies(
                name="euphotic",
                activation_range=(4, 10),
                viability_range=(1, 10),
                production=Production.EXAMPLE[:euphotic],
                diffusion_coefficient=10.0u"m/yr",
                minimum_production=0.01u"m/Myr"),
            M.Facies(
                name="oligophotic",
                production=BenthicProduction(
                    maximum_growth_rate=200.0u"m/Myr",
                    extinction_coefficient=0.1u"m^-1",
                    saturation_intensity=60u"W/m^2"
                ),
                diffusion_coefficient=5.0u"m/yr",
                minimum_production=5.0u"m/Myr"),
            M.Facies(
                name="pelagic",
                active=false,
                production=PelagicProduction(
                    maximum_growth_rate=1.0u"1/Myr",
                    extinction_coefficient=0.1u"m^-1",
                    saturation_intensity=60u"W/m^2"
                ),
                diffusion_coefficient=20.0u"m/yr",
                # minimum_production=10.0u"m/Myr"
            )
        ]

        box = CarboKitten.Box{Periodic{2}}(grid_size=(res, res), phys_scale=phys_scale)

        time_param = TimeProperties(Δt=1.0u"Myr"/steps, steps=steps)

        sea_level(t) =
            10.0u"m" * sin(2π * t / 123456.0u"yr") +
             5.0u"m" * sin(2π * t /  90456.0u"yr")

        input = M.Input(
            time = time_param,
            box = box,
            facies = facies,
            output = output,

            sea_level = sea_level,
            initial_topography = initial_topography,

            insolation = 400.0u"W/m^2",
            subsidence_rate = 30.0u"m/Myr",
            disintegration_rate = 20.0u"m/Myr",
            lithification_time = 100.0u"yr",

            sediment_buffer_size=50,
            depositional_resolution=0.5u"m",

            # diagnostics=true
        )

        run_model(Model{M}, input, "data/output/ca-feedback.h5")
    end
end

Script.main()
# ~/~ end
