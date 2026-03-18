module DenudationParamConfig
using CarboKitten
using CarboKitten.Models: WithDenudation as WDn
using CarboKitten.Production
using Unitful
export facies, input
const m = u"m"
const Myr = u"Myr"
const PERIOD = 0.2Myr
const AMPLITUDE = 20.0m

facies(reactive_surface,mass_density,infiltration_coefficient) = [
    WDn.Facies(
        production=Production.EXAMPLE[:euphotic],
        name="euphotic",
        diffusion_coefficient=50.0u"m/yr",
        reactive_surface=reactive_surface[1],
        mass_density=mass_density[1],
        infiltration_coefficient=infiltration_coefficient[1],
        erodibility = 0.0023u"m/yr"
        ),
    WDn.Facies(
        production=Production.EXAMPLE[:oligophotic],
        name="oligophotic",
        diffusion_coefficient=30.0u"m/yr",
        reactive_surface=reactive_surface[2],
        mass_density=mass_density[2],
        infiltration_coefficient=infiltration_coefficient[2],
        erodibility = 0.0023u"m/yr"
        ),
    WDn.Facies(
        production=Production.EXAMPLE[:aphotic],
        name="aphotic",
        diffusion_coefficient=10.0u"m/yr",
        reactive_surface=reactive_surface[3],
        mass_density=mass_density[3],
        infiltration_coefficient=infiltration_coefficient[3],
        erodibility = 0.0023u"m/yr"
        )
]

facies(erodibility) = [
    WDn.Facies(
        production=Production.EXAMPLE[:euphotic],
        name="euphotic",
        diffusion_coefficient=50.0u"m/yr",
        reactive_surface=10u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = erodibility[1]
        ),
    WDn.Facies(
        production=Production.EXAMPLE[:oligophotic],
        name="oligophotic",
        diffusion_coefficient=30.0u"m/yr",
        reactive_surface=10u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = erodibility[2]
        ),
    WDn.Facies(
        production=Production.EXAMPLE[:aphotic],
        name="aphotic",
        diffusion_coefficient=10.0u"m/yr",
        reactive_surface=10u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = erodibility[3]
        )
]

facies() = [
    WDn.Facies(
        production=Production.EXAMPLE[:euphotic],
        name="euphotic",
        diffusion_coefficient=50.0u"m/yr",
        reactive_surface=10u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.0023u"m/yr"
        ),
    WDn.Facies(
        production=Production.EXAMPLE[:oligophotic],
        name="oligophotic",
        diffusion_coefficient=30.0u"m/yr",
        reactive_surface=10u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.0023u"m/yr"
        ),
    WDn.Facies(
        production=Production.EXAMPLE[:aphotic],
        name="aphotic",
        diffusion_coefficient=10.0u"m/yr",
        reactive_surface=10u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.0023u"m/yr"
        )
]


input(tag, denudation, facies) = WDn.Input(
    tag=tag,
    box=Box{Coast}(grid_size=(100, 50), phys_scale=150.0m),
    time=TimeProperties(
        Δt=0.0002Myr,
        steps=5000),
    output=Dict(
        :topography => OutputSpec(slice=(:,:), write_interval=10),
        :profile => OutputSpec(slice=(:, 25), write_interval=1)),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=t -> AMPLITUDE * sin(2π * t / PERIOD),
    subsidence_rate=20.0m / Myr,
    disintegration_rate=500.0m / Myr,
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5m,
    lithification_time=100.0u"yr",
    facies=facies,
    denudation = denudation)

end