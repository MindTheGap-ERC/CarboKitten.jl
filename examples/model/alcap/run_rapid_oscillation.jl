# Stress test: Rapid oscillation scenario

module ScriptRapidOscillation

using Unitful
using CarboKitten
using CarboKitten.Export: read_slice, data_export, CSV

const PATH = "data/output"
const TAG = "alcap-rapid-oscillation"

const FACIES = [
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=50.0u"m/yr",
        name="euphotic"),
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=25.0u"m/yr",
        name="oligophotic"),
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=12.5u"m/yr",
        name="aphotic"),

    ALCAP.Facies(
        active=false,
        diffusion_coefficient=50.0u"m/yr",
        name="euphotic transported"),
    ALCAP.Facies(
        active=false,
        diffusion_coefficient=50.0u"m/yr",
        name="oligophotic transported"),
    ALCAP.Facies(
        active=false,
        diffusion_coefficient=50.0u"m/yr",
        name="aphotic transported")
]

const PERIOD = 0.1u"Myr"
const AMPLITUDE = 4.5u"m"

const INPUT = ALCAP.Input(
    tag="$TAG",
    box=Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m"),
    time=TimeProperties(
        Δt=0.0002u"Myr",
        steps=5000),
    output=Dict(
        :topography => OutputSpec(slice=(:,:), write_interval=10),
        :profile => OutputSpec(slice=(:, 25), write_interval=1)),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=t -> AMPLITUDE * sin(2π * t / PERIOD),
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=65.0u"m/Myr",  
    disintegration_transfer=p->[0.0u"m", 0.0u"m", 0.0u"m", p[1]+p[4], p[2]+p[5], p[3]+p[6]],
    insolation=400.0u"W/m^2",
    sediment_buffer_size=150,  
    depositional_resolution=1.0u"m",  
    facies=FACIES)

function main()
    run_model(Model{ALCAP}, INPUT, "$(PATH)/$(TAG).h5")
    header, profile = read_slice("$(PATH)/$(TAG).h5", :profile)
    columns = [profile[i] for i in 10:20:70]
    data_export(
        CSV(:sediment_accumulation_curve => "$(PATH)/$(TAG)_sac.csv",
            :age_depth_model => "$(PATH)/$(TAG)_adm.csv",
            :stratigraphic_column => "$(PATH)/$(TAG)_sc.csv",
            :water_depth => "$(PATH)/$(TAG)_wd.csv",
            :metadata => "$(PATH)/$(TAG).toml"),
         header,
         columns)
end

end

ScriptRapidOscillation.main()
