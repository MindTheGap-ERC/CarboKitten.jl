# Model with CA, Production and Active Layer transport (ALCAP)

The following **S**edimentation model includes the Burgess 2013 **C**ellular **A**utomaton, Bosscher & Schlager 1992 **P**roduction curves and an **A**ctive **L**ayer transport model, based on Paola 1992, henceforth ALCAP.

![Result from alternative input](fig/alcaps-alternative.png)

## Example Input

The following is a complete example input.

``` {.julia #alcap-example-input}
const TAG = "alcap-example"

const FACIES = [
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=50.0u"m/yr"),
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=25.0u"m/yr"),
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=12.5u"m/yr")
]

const PERIOD = 0.2u"Myr"
const AMPLITUDE = 4.0u"m"

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
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=0.5u"m",
    facies=FACIES)
```

``` {.julia .task file=examples/model/alcap/run.jl}
#| requires: src/Models/ALCAP.jl
#| creates: data/output/alcap-example.h5

module Script

using Unitful
using CarboKitten
using CarboKitten.Export: read_slice, data_export, CSV

const PATH = "data/output"

<<alcap-example-input>>

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

Script.main()
```

```@raw html
<details><summary>Plotting code</summary>
```

``` {.julia .task file=examples/model/alcap/plot.jl}
#| creates: docs/src/_fig/alcaps-alternative.png
#| requires: data/output/alcap-example.h5
#| collect: figures

using GLMakie
using CarboKitten.Visualization

GLMakie.activate!()

save("docs/src/_fig/alcaps-alternative.png", summary_plot("data/output/alcap-example.h5"))
```

```@raw html
</details>
```

## Modular Implementation

```component-dag
CarboKitten.Models.ALCAP
```

``` {.julia file=src/Models/ALCAP/Example.jl}
module Example

using Unitful
using ..ALCAP: ALCAP
using CarboKitten.Boxes: Box, Coast
using CarboKitten.Config: TimeProperties
using CarboKitten.OutputData: OutputSpec

<<alcap-example-input>>

end
```

``` {.julia file=src/Models/ALCAP.jl}
@compose module ALCAP
@mixin Tag, H5Writer, CAProduction, ActiveLayer, InitialSediment

using ..Common
using ..CAProduction: production
using ..TimeIntegration
using ..WaterDepth: water_depth
using ModuleMixins: @for_each

export Input, Facies

function initial_state(input::AbstractInput)
    ca_state = CellularAutomaton.initial_state(input)
    for _ in 1:20
        CellularAutomaton.step!(input)(ca_state)
    end

    sediment_height = zeros(Height, input.box.grid_size...)
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...)
    active_layer = zeros(Amount, n_facies(input), input.box.grid_size...)

    state = State(
        step=0, sediment_height=sediment_height,
        sediment_buffer=sediment_buffer,
        active_layer = active_layer,
        ca=ca_state.ca, ca_priority=ca_state.ca_priority)

    InitialSediment.push_initial_sediment!(input, state)
    return state
end

function step!(input::Input)
    step_ca! = CellularAutomaton.step!(input)
    disintegrate! = ActiveLayer.disintegrator(input)
    produce = production(input)
    transport! = ActiveLayer.transporter(input)
    local_water_depth = water_depth(input)
    na = [CartesianIndex()]
    pf = cementation_factor(input)

    function (state::State)
        if mod(state.step, input.ca_interval) == 0
            step_ca!(state)
        end

        wd = local_water_depth(state)
        p = produce(state, wd)
        d = disintegrate!(state)

        state.active_layer .+= p
        state.active_layer .+= d
        transport!(state)

        deposit = pf .* state.active_layer
        push_sediment!(state.sediment_buffer, deposit ./ input.depositional_resolution .|> NoUnits)
        state.active_layer .-= deposit
        state.sediment_height .+= sum(deposit; dims=1)[1,:,:]
        state.step += 1

        return Frame(
            production = p,
            disintegration = d,
            deposition = deposit)
    end
end

function write_header(fid, input::AbstractInput)
    @for_each(P -> P.write_header(fid, input), PARENTS)
end

include("ALCAP/Example.jl")

end
```

## API

```@autodocs
Modules = [CarboKitten.Models.ALCAP]
```
