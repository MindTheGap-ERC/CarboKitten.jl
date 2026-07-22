# ~/~ begin <<docs/src/models/without-ca.md#src/Models/WithoutCA.jl>>[init]
@compose module WithoutCA
@mixin Tag, Diagnostics, Output, Production, ActiveLayer, InitialSediment

using ..Common
using ..Production: uniform_production
using ..TimeIntegration
using ..WaterDepth
using ...Output: Frame
using ModuleMixins: @for_each

export Input, Facies

function initial_state(input::Input)
    sediment_thickness = zeros(Height, input.box.grid_size...)
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...)
    active_layer = zeros(Amount, n_facies(input), input.box.grid_size...)
    state = State(
        step=0,
        bathymetry=initial_topography(input),
        sediment_thickness=sediment_thickness,
        sediment_buffer=sediment_buffer,
        active_layer=active_layer)
    InitialSediment.push_initial_sediment!(input, state)
    return state
end

function initial_frame(input::Input)
    dep = stack(InitialSediment.initial_sediment(input.box, f) for f in input.facies; dims=1)
    return Frame(production=zeros(Sediment,size(dep)),
                  disintegration=zeros(Sediment,size(dep)),
                  deposition=dep)
end

function step!(input::Input)
    disintegrate! = ActiveLayer.disintegrator(input)
    transport! = ActiveLayer.transporter(input)
    produce = uniform_production(input)
    dt = input.time.Δt
    local_water_depth = water_depth(input)
    na = [CartesianIndex()]
    pf = lithification_factor(input)
    dtf = input.disintegration_transfer
    push! = push_sediment(input)
    subside! = subsider(input)

    function (state::State)
        wd = local_water_depth(state)
        p = produce(state, wd)
        d = disintegrate!(state)

        state.active_layer .+= p
        state.active_layer .+= dtf(d)
        transport!(state)

        deposit = pf .* state.active_layer
        push!(state, deposit)
        state.active_layer .-= deposit
        subside!(state)
        state.step += 1

        return Frame(
            production=p,
            disintegration=d,
            deposition=deposit)
    end
end

function write_header(input::AbstractInput, output::AbstractOutput)
    @for_each(P -> P.write_header(input, output), PARENTS)
end

end
# ~/~ end
