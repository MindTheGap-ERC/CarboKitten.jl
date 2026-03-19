# ~/~ begin <<docs/src/components/cellular-automata.md#examples/model/alcap/ca-feedback.jl>>[init]
module Script
    using CarboKitten
    using CarboKitten.Production
    using CarboKitten.Models: ALCAP as M

    initial_topography(x, y) =
        min(0.0u"m", - sqrt((x - 7.5u"km")^2 + (y - 7.5u"km")^2) / 100.0 + 20.0u"m")

    function main()
        res = 100
        steps = 5000
        phys_scale = 15.0u"km" / res

        output = Dict(
            :topography => OutputSpec(write_interval = max(1, div(steps, 50))),
            :profile    => OutputSpec(slice = (:, div(res, 2)+1)))

        facies(feedback) = [
            M.Facies(
                name="euphotic",
                activation_range=(4, 10),
                viability_range=(1, 10),
                production=Production.EXAMPLE[:euphotic],
                transport_coefficient=10.0u"m/yr",
                minimum_production=feedback ? 0.01u"m/Myr" : nothing),
            M.Facies(
                name="oligophotic",
                production=BenthicProduction(
                    maximum_growth_rate=200.0u"m/Myr",
                    extinction_coefficient=0.1u"m^-1",
                    saturation_intensity=60u"W/m^2"
                ),
                transport_coefficient=5.0u"m/yr",
                minimum_production=feedback ? 5.0u"m/Myr" : nothing),
            M.Facies(
                name="pelagic",
                active=false,
                production=PelagicProduction(
                    maximum_growth_rate=1.0u"1/Myr",
                    extinction_coefficient=0.1u"m^-1",
                    saturation_intensity=60u"W/m^2"
                ),
                transport_coefficient=20.0u"m/yr",
                # minimum_production=10.0u"m/Myr"
            )
        ]

        box = CarboKitten.Box{Periodic{2}}(grid_size=(res, res), phys_scale=phys_scale)

        time_param = TimeProperties(Δt=1.0u"Myr"/steps, steps=steps)

        sea_level(t) =
            10.0u"m" * sin(2π * t / 123456.0u"yr") +
             5.0u"m" * sin(2π * t /  80456.0u"yr")

        input(feedback) = M.Input(
            time = time_param,
            box = box,
            facies = facies(feedback),
            output = output,

            sea_level = sea_level,
            initial_topography = initial_topography,
            ca_interval = 10,

            insolation = 400.0u"W/m^2",
            subsidence_rate = 30.0u"m/Myr",
            disintegration_rate = 20.0u"m/Myr",
            lithification_time = 100.0u"yr",

            sediment_buffer_size=50,
            depositional_resolution=0.5u"m",

            # diagnostics=true
        )

        run_model(Model{M}, input(false), "data/output/ca-wo-feedback.h5")
        run_model(Model{M}, input(true), "data/output/ca-feedback.h5")
    end
end

Script.main()
# ~/~ end
