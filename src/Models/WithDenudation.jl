# ~/~ begin <<docs/src/models/with-denudation.md#src/Models/WithDenudation.jl>>[init]
@compose module WithDenudation
@mixin Tag, Diagnostics, Output, CAProduction, ActiveLayer, Denudation

using ..Common
using ..CAProduction: production
using ..TimeIntegration
using ..WaterDepth
using ModuleMixins: @for_each
using ...Stencil
using ...BoundaryTrait
using ...Denudation.EmpiricalDenudationMod: slope_kernel
using ...Output: Frame
export Input, Facies

function initial_state(input::AbstractInput)
    ca_state = CellularAutomaton.initial_state(input)
    for _ in 1:20
        CellularAutomaton.step!(input)(ca_state)
    end

    sediment_height = zeros(Height, input.box.grid_size...)
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...)
    active_layer = zeros(Amount, n_facies(input), input.box.grid_size...)

    return State(
        step=0, sediment_height=sediment_height,
        sediment_buffer=sediment_buffer,
         active_layer=active_layer,
        ca=ca_state.ca, ca_priority=ca_state.ca_priority)
end

function step!(input::Input)
    step_ca! = CellularAutomaton.step!(input)
    disintegrate! = ActiveLayer.disintegrator(input)
    produce = production(input)
    transport! = ActiveLayer.transporter(input)
    denudate = denudation(input)
    redistribute = redistribution(input)
    local_water_depth = water_depth(input)
    slopefn = slope_function(input, input.box)
    pf = cementation_factor(input)
    dtf = input.disintegration_transfer

    slope = Array{Float64}(undef, input.box.grid_size...)
    denuded_sediment = Array{Float64}(undef, n_facies(input), input.box.grid_size...)

    function (state::State)
        if mod(state.step, input.ca_interval) == 0
            step_ca!(state)
        end

        wd = local_water_depth(state)
        w = wd ./ u"m"
        slopefn(w, slope, input.box.phys_scale ./ u"m")

        p = produce(state, wd)
        d = disintegrate!(state)

        state.active_layer .+= p
        state.active_layer .+= dtf(d)

        transport!(state)

        deposit = pf .* state.active_layer
        push_sediment!(state.sediment_buffer, deposit ./ input.depositional_resolution .|> NoUnits)
        state.active_layer .-= deposit
        state.sediment_height .+= sum(deposit; dims=1)[1, :, :]

        # subaerial: denudation and redistribution
        denudation_mass = denudate(state, w, slope)
        if denudation_mass !== nothing
            denudation_mass = denudation_mass |> x -> sum(x, dims=1) |> x -> dropdims(x, dims=1) |> x -> min.(x, state.sediment_height)

            state.sediment_height .-= denudation_mass
            pop_sediment!(state.sediment_buffer, denudation_mass ./ input.depositional_resolution .|> NoUnits, denuded_sediment)

            d .+= denuded_sediment .* input.depositional_resolution

            redistribution_mass = redistribute(state, w, denuded_sediment .* input.depositional_resolution)
            if redistribution_mass !== nothing
                # redistribution returns a 3D facies array in meters
                push_sediment!(state.sediment_buffer, redistribution_mass ./ input.depositional_resolution .|> NoUnits)
                state.sediment_height .+= sum(redistribution_mass; dims=1)[1, :, :]
            end
        end

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
