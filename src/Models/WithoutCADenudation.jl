# ~/~ begin <<docs/src/models/without-ca.md#src/Models/WithoutCA.jl>>[init]
module WithoutCADenudation
import ModuleMixins.Procedural: requires
import ..Tag, ..Diagnostics, ..Output, ..Production, ..ActiveLayer, ..Denudation, ..InitialSediment

using ..Common
using ..FaciesBase: n_facies
using ..SedimentBuffer: push_sediment, pop_sediment
using ..ActiveLayer: lithification_factor
using ..Denudation: slope_function, denudation, redistribution
using ..Production: uniform_production
using ..TimeIntegration
using ..WaterDepth
using ...Stencil
using ...BoundaryTrait
using ...Output: Frame
using ...Denudation.EmpiricalDenudationMod: slope_kernel

# The model composes the records of its components. These are defined in the
# model module itself, so `requires` is defined first, instead of using
# `@component`, which defines it at the end of the module.
requires(::Val{WithoutCADenudation}) = (Tag, Diagnostics, Output, Production, ActiveLayer, Denudation, InitialSediment)
const Input = record_type(:Input, WithoutCADenudation)
const State = record_type(:State, WithoutCADenudation)
const Facies = record_type(:Facies, WithoutCADenudation)

export Input, Facies

initial_state(input::Input) = Common.initial_state(State, input)

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
    local_water_depth = water_depth(input)
    slopefn = slope_function(input, input.box)
    pf = lithification_factor(input)
    dtf = input.disintegration_transfer
    push! = push_sediment(input)
    pop! = pop_sediment(input)
    subside! = subsider(input)
    
    slope = Array{Float64}(undef, input.box.grid_size...)
    denuded_sediment = Array{Amount, 3}(undef, n_facies(input), input.box.grid_size...)

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
        push!(state, deposit) 
        state.active_layer .-= deposit
        
        # denudation and redistribution 
        denudation_mass = denudate(state, w, slope)
        if denudation_mass !== nothing
            denudation_mass = denudation_mass |>
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

end
# ~/~ end
