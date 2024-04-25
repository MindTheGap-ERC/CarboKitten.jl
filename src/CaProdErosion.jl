# ~/~ begin <<docs/src/ca-prod-with-erosion.md#src/CaProdErosion.jl>>[init]
module CaProdErosion

using CarboKitten
using ..Stencil: Periodic, stencil
using ..Utility
#using ..BS92: sealevel_curve
using ..Denudation: DenudationType, denudation, DenudationFrame
using ..Burgess2013

using HDF5
using .Iterators: drop, peel, partition, map, take

# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-input>>[init]
@kwdef struct Facies
    viability_range::Tuple{Int, Int}
    activation_range::Tuple{Int, Int}

    maximum_growth_rate::typeof(1.0u"m/Myr")
    extinction_coefficient::typeof(1.0u"m^-1")
    saturation_intensity::typeof(1.0u"W/m^2")

    reactive_surface::Float64 #reactive surface
    mass_density::Float64 #density of different carb factory
    infiltration_coefficient::Float64 #infiltration coeff
end

@kwdef struct Input
    box :: Box
    time :: TimeProperties

    sea_level       # Myr -> m
    subsidence_rate::typeof(1.0u"m/Myr")
    initial_depth   # m -> m

    facies::Vector{Facies}
    insolation::typeof(1.0u"W/m^2")

    denudation::DenudationType
end
# ~/~ end
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-frame>>[init]

abstract type Frame end
struct ProductionFrame <: Frame
    production::Array{Float64,3}
end

struct OutputFrame
    production::Array{Float64,3}
    denudation::Array{Float64,2}
    redistribution::Array{Float64,2}
end
# ~/~ end
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-state>>[init]
# FIXME: deduplicate
mutable struct State
    time::Float64
    height::Array{Float64,2}
end
# ~/~ end
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-model>>[init]
# FIXME: deduplicate
function initial_state(input::Input)  # -> State
    height = zeros(Float64, input.grid_size...)
    for i in CartesianIndices(height)
        height[i] = input.initial_depth(i[2] * input.phys_scale)
    end
    return State(0.0, height)
end
# ~/~ end
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-model>>[1]
# FIXME: deduplicate
function propagator(input::Input)
    # ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-init-propagator>>[init]
    n_facies = length(input.facies)
    ca_init = rand(0:n_facies, input.grid_size...)
    ca = drop(run_ca(Periodic{2}, input.facies, ca_init, 3), 20)

    function water_depth(s::State)
        s.height .- input.sea_level(s.time)
    end
    # prepare functions for erosion
    # ~/~ end
    slopefn = stencil(Float64, Periodic{2}, (3, 3), slope_kernel)
    function (s::State)  # -> Frame
        # ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-propagate>>[init]
        production = zeros(Float64, input.grid_size..., n_facies)
        denudation = zeros(Float64, input.grid_size...)
        redistribution = zeros(Float64, input.grid_size...)
        slope = zeros(Float64, input.grid_size...)
        facies_map, ca = peel(ca)
        w = water_depth(s)
        slopefn(w,slope,input.phys_scale) # slope is calculated with square so no need for -w
        if input.erosion_type == 2
        redis = mass_erosion(Float64,Periodic{2},slope,(3,3),w,input.phys_scale,input.facies.inf)
        redistribution = total_mass_redistribution(redis, slope)
        else
        redistribution = redistribution
        end
        Threads.@threads for idx in CartesianIndices(facies_map)
            f = facies_map[idx]
                if f == 0
                    continue
                end
            if w[idx] > 0.0
                production[Tuple(idx)..., f] = production_rate(input.insolation, input.facies[f], w[idx])
            else
                if input.erosion_type == 1
                    denudation[Tuple(idx)...] = dissolution(input.temp,input.precip,input.alpha,input.pco2,w[idx],input.facies[f])
                elseif input.erosion_type == 2
                    denudation[Tuple(idx)...] = physical_erosion(slope[idx],input.facies.inf)
                elseif input.erosion_type == 3
                    denudation[Tuple(idx)...] = emperical_denudation(input.precip, slope[idx])
                elseif input.erosion_type == 0
                    denudation[Tuple(idx)...] = denudation[Tuple(idx)...]
                end
            end
        end
        return Frame(production, denudation, redistribution)#
        # ~/~ end
    end
end
# ~/~ end
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-model>>[2]
function updater(input)
    n_facies = length(input.facies)
    function update(state, Δ::ProductionFrame)
        state.height .-= sum(Δ.production; dims=3) .* input.Δt
        state.height .+= Δ.denudation .* input.Δt  #number already in kyr
        state.height .-= Δ.redistribution .* input.Δt
        state.height .+= input.subsidence_rate * input.Δt
        state.time += input.Δt
    end

    function update(state, Δ::DenudationFrame)
        # FIXME: implement
    end

    update
end
# ~/~ end
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-model>>[3]
function run_model(input::Input)
    Channel{OutputFrame}() do ch
        s = initial_state(input)
        p = Production::propagator(input)
        d = Denudation::propagator(input)  # FIXME: implement
        u = updater(input)

        while true
            Δ_prod = p(s)
            u(s, Δ_prod)
            Δ_denu = d(s)
            u(s, Δ_denu)
            put!(ch, OutputFrame(Δ_prod.production, Δ_denu.denudation, Δ_denu.redistribution))
        end
    end
end
# ~/~ end

function stack_frames(fs::Vector{OutputFrame})  # -> Frame
    OutputFrame(sum(f.production for f in fs),sum(f.denudation for f in fs),sum(f.redistribution for f in fs))#
end

function main(input::Input, output::String)
    x_axis = (0:(input.grid_size[2]-1)) .* input.phys_scale
    y_axis = (0:(input.grid_size[1]-1)) .* input.phys_scale
    initial_height = input.initial_depth.(x_axis)
    n_writes = input.time_steps ÷ input.write_interval

    h5open(output, "w") do fid
        gid = create_group(fid, "input")
        gid["x"] = collect(x_axis)
        gid["y"] = collect(y_axis)
        gid["height"] = collect(initial_height)
        gid["t"] = collect((0:(n_writes-1)) .* (input.Δt * input.write_interval))
        attr = attributes(gid)
        attr["delta_t"] = input.Δt
        attr["write_interval"] = input.write_interval
        attr["time_steps"] = input.time_steps
        attr["subsidence_rate"] = input.subsidence_rate

        n_facies = length(input.facies)
        ds = create_dataset(fid, "sediment", datatype(Float64),
            dataspace(input.grid_size..., n_facies, input.time_steps),
            chunk=(input.grid_size..., n_facies, 1))
        denudation = create_dataset(fid, "denudation", datatype(Float64),
            dataspace(input.grid_size..., input.time_steps),
           chunk=(input.grid_size..., 1))
        redistribution = create_dataset(fid, "redistribution", datatype(Float64),
           dataspace(input.grid_size..., input.time_steps),
          chunk=(input.grid_size..., 1))

        results = map(stack_frames, partition(run_model(input), input.write_interval))
        for (step, frame) in enumerate(take(results, n_writes))
            ds[:, :, :, step] = frame.production
            denudation[:,:,step] = frame.denudation
            redistribution[:,:,step] = frame.redistribution
        end
    end
end

end # CaProd
# ~/~ end