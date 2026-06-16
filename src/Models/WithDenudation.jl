# ~/~ begin <<docs/src/models/with-denudation.md#src/Models/WithDenudation.jl>>[init]
@compose module WithDenudation
@mixin Tag, Diagnostics, Output, CAProduction, ActiveLayer, Denudation, InitialSediment

using ..Common
using ..CAProduction: production
using ..TimeIntegration
using ..WaterDepth
using ModuleMixins: @for_each
using ...Stencil
using ...BoundaryTrait
using ...Denudation.EmpiricalDenudationMod: slope_kernel
using ...Output: Frame
export Input, Facies, BenthicProduction, PelagicProduction

function initial_state(input::AbstractInput)
    ca_state = CellularAutomaton.initial_state(input)
    for _ in 1:20
        CellularAutomaton.step!(input)(ca_state)
    end

    bathymetry = initial_topography(input)
    sediment_thickness = zeros(Height, input.box.grid_size...)
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...)
    active_layer = zeros(Amount, n_facies(input), input.box.grid_size...)

    state = State(
        step=0,
        bathymetry=bathymetry,
        sediment_thickness=sediment_thickness,
        sediment_buffer=sediment_buffer,
        active_layer=active_layer,
        ca=ca_state.ca, ca_priority=ca_state.ca_priority)

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
    step_ca! = CellularAutomaton.step!(input)
    disintegrate! = ActiveLayer.disintegrator(input)
    produce = production(input)
    transport! = ActiveLayer.transporter(input)
    denudate = denudation(input)
    redistribute = redistribution(input)
    local_water_depth = water_depth(input)
    slopefn = slope_function(input, input.box)
    pf = lithification_factor(input)
    dtf = input.disintegration_transfer
    pop! = pop_sediment!(input)
    push! = push_sediment(input)
    subside! = subsider(input)

    slope = Array{Float64}(undef, input.box.grid_size...)
    denuded_sediment = Array{Amount}(undef, n_facies(input), input.box.grid_size...)

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
        push!(state, deposit)
        state.active_layer .-= deposit

        # subaerial: denudation and redistribution
        denudation_mass = denudate(state, w, slope)
        if denudation_mass !== nothing
            denudation_mass = denudation_mass |>
                x -> sum(x, dims=1) |>
                x -> dropdims(x, dims=1) |>
                x -> min.(x, state.sediment_thickness)
            pop!(state, denudation_mass, denuded_sediment)

            d .+= denuded_sediment

            redistribution_mass = redistribute(state, w, denuded_sediment)
            if redistribution_mass !== nothing
                # redistribution returns a 3D facies array in meters
                push!(state, redistribution_mass)
            end
        end

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
