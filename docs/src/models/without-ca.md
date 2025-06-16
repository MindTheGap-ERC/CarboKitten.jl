# Model without CA

The following model is similar to the ALCAP model, with the exception that production is not managed by the cellular automaton, rather each facies produces in all grid cells according to their own production curves.

``` {.julia file=src/Models/WithoutCA.jl}
@compose module WithoutCA
@mixin Tag, H5Writer, Production, ActiveLayer

using ..Common
using ..Production: uniform_production
using ..TimeIntegration
using ..WaterDepth
using ModuleMixins: @for_each

export Input, Facies

function initial_state(input::Input)
    sediment_height = zeros(Height, input.box.grid_size...)
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...)
    return State(step=0, sediment_height=sediment_height, sediment_buffer=sediment_buffer)
end

function step!(input::Input)
    disintegrate! = ActiveLayer.disintegrator(input)
    transport! = ActiveLayer.transporter(input)
    τ = uniform_production(input)
    dt = input.time.Δt

    function (state::State)
        p = τ(state) .* dt
        d = disintegrate!(state)

        active_layer = p .+ d
        transport!(state, active_layer)

        push_sediment!(state.sediment_buffer, active_layer ./ input.depositional_resolution .|> NoUnits)
        state.sediment_height .+= sum(active_layer; dims=1)[1,:,:]
        state.step += 1

        return Frame(
            production = p,
            disintegration = d,
            deposition = active_layer)
    end
end

function write_header(fid, input::AbstractInput)
    @for_each(P -> P.write_header(fid, input), PARENTS)
end

end
```
