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
                diffusion_coefficient=25.0u"m/yr"),
            M.Facies(
                name="oligophotic",
                production=Production.EXAMPLE[:oligophotic],
                diffusion_coefficient=10.0u"m/yr"),
            M.Facies(
                name="pelagic",
                active=false,
                production=Production.EXAMPLE[:pelagic],
                diffusion_coefficient=50.0u"m/yr")
        ]
        box = CarboKitten.Box(Periodic{2}, grid_size=(256, 256), phys_scale=res)
        time_param = TimeParameters(Δt=0.0002u"Myr", steps=steps),
        input = M.Input(
            time = time_param,
            box = box,
            facies = facies,

            sea_level=t -> AMPLITUDE * sin(2π * t / PERIOD),
            initial_topography = initial_topography,

            insolation = 400.0u"W/m^2",
            subsidence_rate = 50.0u"m/Myr",
            disintegration_rate = 50.0u"m/Myr",
            lithification_time = 100.0u"yr",

            sediment_buffer_size=50,
            depositional_resolution=0.5u"m",
        )
    end
end

Script.main()
# ~/~ end
