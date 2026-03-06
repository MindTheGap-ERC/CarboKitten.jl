# ~/~ begin <<docs/src/components/cellular-automata.md#examples/model/alcap/ca-feedback.jl>>[init]
module Script
    using CarboKitten
    using CarboKitten.Production
    using CarboKitten.Models: ALCAP as M

    initial_topography(x, y) =
        - sqrt((x - 7.5u"km")^2 + (y - 7.5u"km")^2) / 100.0

    function main()
        res = 150.0u"m"
        steps = 5000

        output(res, steps) = Dict(
            :topography => OutputSpec(write_interval = max(1, div(steps, 50))),
            :profile    => OutputSpec(slice = (:, div(res, 2)+1)))

        facies = [
            M.Facies(
                name="euphotic",
                production=Production.EXAMPLE[:euphotic],
                diffusion_coefficient=5.0u"m/yr"),
            M.Facies(
                name="oligophotic",
                production=Production.EXAMPLE[:oligophotic],
                diffusion_coefficient=1.0u"m/yr"),
            M.Facies(
                name="aphotic",
                # active=false,
                production=Production.EXAMPLE[:aphotic],
                diffusion_coefficient=2.0u"m/yr")
        ]

        box = CarboKitten.Box{Periodic{2}}(grid_size=(256, 256), phys_scale=res)

        time_param = TimeProperties(Δt=1.0u"Myr"/steps, steps=steps)

        sea_level(t) =
            10.0u"m" * sin(2π * t / 123456.0u"yr") +
             5.0u"m" * sin(2π * t /  23456.0u"yr")

        input = M.Input(
            time = time_param,
            box = box,
            facies = facies,

            sea_level = sea_level,
            initial_topography = initial_topography,

            insolation = 400.0u"W/m^2",
            subsidence_rate = 50.0u"m/Myr",
            disintegration_rate = 50.0u"m/Myr",
            lithification_time = 100.0u"yr",

            sediment_buffer_size=50,
            depositional_resolution=0.5u"m",
        )

        run_model(Model{M}, input, "data/output/ca-feedback.h5")
    end
end

Script.main()
# ~/~ end
