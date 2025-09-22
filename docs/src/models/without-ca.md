# Model without CA

The following model is similar to the ALCAP model, with the exception that production is not managed by the cellular automaton, rather each facies produces in all grid cells according to their own production curves.

``` {.julia file=src/Models/WithoutCA.jl}
@compose module WithoutCA
@mixin Tag, Output, Production, ActiveLayer

using ..Common
using ..Production: uniform_production
using ..TimeIntegration
using ..WaterDepth
using ModuleMixins: @for_each

export Input, Facies

function initial_state(input::Input)
    sediment_height = zeros(Height, input.box.grid_size...)
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...)
    active_layer = zeros(Amount, n_facies(input), input.box.grid_size...)
    return State(step=0, sediment_height=sediment_height, sediment_buffer=sediment_buffer, active_layer=active_layer)
end

function step!(input::Input)
    disintegrate! = ActiveLayer.disintegrator(input)
    transport! = ActiveLayer.transporter(input)
    produce = uniform_production(input)
    dt = input.time.Δt
    local_water_depth = water_depth(input)
    na = [CartesianIndex()]
    pf = cementation_factor(input)

    function (state::State)
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

end
```
