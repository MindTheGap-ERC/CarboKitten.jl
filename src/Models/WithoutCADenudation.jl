# ~/~ begin <<docs/src/models/without-ca.md#src/Models/WithoutCA.jl>>[init]
@compose module WithoutCADenudation
@mixin Tag, Diagnostics, Output, Production, ActiveLayer, Denudation, InitialSediment

using ..Common
using ..Production: uniform_production
using ..TimeIntegration
using ..WaterDepth
using ...Output: Frame
using ModuleMixins: @for_each
using ...Denudation.EmpiricalDenudationMod: slope_kernel

export Input, Facies

function initial_state(input::Input)
    sediment_height = zeros(Height, input.box.grid_size...)
    sediment_buffer = zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...)
    active_layer = zeros(Amount, n_facies(input), input.box.grid_size...)
    state = State(step=0, sediment_height=sediment_height, sediment_buffer=sediment_buffer, active_layer=active_layer)
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
    denudate = denudation(input)
    redistribute = redistribution(input)
    dt = input.time.Δt
    local_water_depth = water_depth(input)
    na = [CartesianIndex()]
    pf = lithification_factor(input)
    dtf = input.disintegration_transfer
    
    slopefn = slope_function(input, input.box)
    slope = Array{Float64}(undef, input.box.grid_size...)
    denuded_sediment = Array{Float64}(undef, n_facies(input), input.box.grid_size...)

    function (state::State)
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
        
        # denudation and redistribution 
        denudation_mass = denudate(state, w, slope)
        if denudation_mass !== nothing
            denudation_mass = denudation_mass |> x -> min.(x, state.sediment_height)

            state.sediment_height .-= denudation_mass

            pop_sediment!(state.sediment_buffer, denudation_mass ./ input.depositional_resolution .|> NoUnits, denuded_sediment)

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
