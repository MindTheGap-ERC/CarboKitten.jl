# Example of a custom model: active-layer transport with a user-defined
# production function `production(t, x, y, water_depth)` that returns a
# production rate for each facies.

module CustomProductionModel

import ModuleMixins.Procedural: requires
using CarboKitten
using CarboKitten.Components: Tag, Diagnostics, Output, ActiveLayer, InitialSediment
using CarboKitten.Components.Common
using CarboKitten.Components.FaciesBase: n_facies
using CarboKitten.Components.SedimentBuffer: push_sediment
using CarboKitten.Components.ActiveLayer: lithification_factor
using CarboKitten.Components.TimeIntegration: TimeIntegration
using CarboKitten.Components.WaterDepth: water_depth, subsider
using CarboKitten.Output: Frame

# A minimal component that adds the `production` field to the input.
module CustomProduction
using CarboKitten.Components.Common

@kwdef struct Input <: InputPart
    production   # function of (t, x, y, water_depth), returning a vector of rates
end
end

requires(::Val{CustomProductionModel}) = (Tag, Diagnostics, Output, ActiveLayer, InitialSediment, CustomProduction)
const Input = record_type(:Input, CustomProductionModel)
const State = record_type(:State, CustomProductionModel)
const Facies = record_type(:Facies, CustomProductionModel)

initial_state(input::Input) = Common.initial_state(input)

function initial_frame(input::Input)
    dep = stack(InitialSediment.initial_sediment(input.box, f) for f in input.facies; dims=1)
    return Frame(production=zeros(Sediment, size(dep)),
                 disintegration=zeros(Sediment, size(dep)),
                 deposition=dep)
end

function step!(input::Input)
    disintegrate! = ActiveLayer.disintegrator(input)
    transport! = ActiveLayer.transporter(input)
    local_water_depth = water_depth(input)
    get_time = TimeIntegration.time(input)
    x, y = box_axes(input.box)
    dt = input.time.Δt
    pf = lithification_factor(input)
    dtf = input.disintegration_transfer
    push! = push_sediment(input)
    subside! = subsider(input)

    p = Array{Amount, 3}(undef, n_facies(input), input.box.grid_size...)
    deposit = Array{Amount, 3}(undef, n_facies(input), input.box.grid_size...)

    function (state::State)
        t = get_time(state)
        wd = local_water_depth(state)
        for j in eachindex(y), i in eachindex(x)
            p[:, i, j] .= input.production(t, x[i], y[j], wd[i, j]) .* dt
        end
        d = disintegrate!(state)

        state.active_layer .+= p
        state.active_layer .+= dtf(d)
        transport!(state)

        deposit .= pf .* state.active_layer
        push!(state, deposit)
        state.active_layer .-= deposit
        subside!(state)
        state.step += 1

        return Frame(
            production=copy(p),
            disintegration=d,
            deposition=copy(deposit))
    end
end

function main()
    # a single facies producing at a fixed rate in water shallower than 20 m
    production(t, x, y, wd) = [0.0 < wd / u"m" < 20.0 ? 200.0u"m/Myr" : 0.0u"m/Myr"]

    input = Input(
        box=Box{Coast}(grid_size=(100, 1), phys_scale=150.0u"m"),
        time=TimeProperties(Δt=200.0u"yr", steps=1000),
        output=Dict(:profile => OutputSpec(slice=(:, 1))),
        facies=[Facies(transport_coefficient=10.0u"m/yr")],
        production=production,
        initial_topography=(x, y) -> -x / 300.0,
        sea_level=t -> 0.0u"m",
        subsidence_rate=50.0u"m/Myr",
        disintegration_rate=50.0u"m/Myr")

    run_model(Model{CustomProductionModel}, input, MemoryOutput(input))
end

end

CustomProductionModel.main()
