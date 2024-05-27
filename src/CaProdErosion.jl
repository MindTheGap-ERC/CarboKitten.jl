# ~/~ begin <<docs/src/ca-prod-with-erosion.md#src/CaProdErosion.jl>>[init]
module CaProdErosion

using CarboKitten
using ..Stencil: Periodic, stencil
using ..Utility
#using ..BS92: sealevel_curve
using ..Denudation: denudation
using ..Burgess2013
using ..InputConfig: Input, DenudationType
using Unitful
using HDF5
using .Iterators: drop, peel, partition, map, take
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-input>>[init]
struct Facies
    viability_range::Tuple{Int, Int}
    activation_range::Tuple{Int, Int}

    maximum_growth_rate::typeof(1.0u"m/Myr")
    extinction_coefficient::typeof(1.0u"m^-1")
    saturation_intensity::typeof(1.0u"W/m^2")

    reactive_surface::Float64 #reactive surface
    mass_density::Float64 #density of different carb factory
    infiltration_coefficient::Float64 #infiltration coeff
end

# ~/~ end
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-frame>>[init]

abstract type Frame end
struct ProductionFrame <: Frame
    production::Array{Float64,3}
end

struct DenudationFrame <: Frame
    denudation::Array{Float64,2}
    redistribution::Array{Float64,2}
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

# propagator for production
function prod_propagator(input::Input)
    n_facies = length(input.facies)
    ca_init = rand(0:n_facies, input.grid_size...)
    ca = drop(run_ca(Periodic{2}, input.facies, ca_init, 3), 20)

    function water_depth(s::State)
        s.height .- input.sea_level(s.time)
    end

    function(s::State)
    production = zeros(typeof(0.0u"m/Myr"), input.box.grid_size..., n_facies)
    facies_map, ca = peel(ca)
    w = water_depth(s)
    for idx in CartesianIndices(facies_map)
        f = facies_map[idx]
        if f == 0
            continue
        end
        production[Tuple(idx)..., f] = production_rate(input.insolation, input.facies[f], w[idx])
    end
    end

    return ProductionFrame(production)
end


function denu_propagator(input::Input)

    function water_depth(state)
        state.height .- state.sea_level
    end

    w = water_depth(s)

    function (s::State)
        denudation = zeros(input.box.grid_size...)
        redistribution = zeros(input.box.grid_size...)
    for idx in CartesianIndices(w)
        if w[idx] < 0.0
            (denudation[Tuple(idx)],redistribution[Tuple(idx)]) = denudation(input,input.denu_param,s)
        end
    end

    end
    return DenudationFrame(denudation,redistribution)
end

# ~/~ end
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-model>>[2]
function updater(input)
    n_facies = length(input.facies)
    function update(state, Δ::ProductionFrame)
        state.height .-= sum(Δ.production; dims=3) .* input.Δt
        state.height .+= Δ.denudation .* input.Δt  #number already in kyr
        state.time += input.Δt
    end

    function update(state, Δ::DenudationFrame)
        # FIXME: implement
        state.hieght .+= Δ.denudation .* input.Δt
        state.height .-= Δ.redistribution .* input.Δt
        state.time += input.Δt
    end

    update
end
# ~/~ end
# ~/~ begin <<docs/src/ca-prod-with-erosion.md#cape-model>>[3]
function run_model(input::Input)
    Channel{OutputFrame}() do ch
        s = initial_state(input)
        p = prod_propagator(input)
        d = denu_propagator(input)  # FIXME: implement
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
            ds[:, :, :, step] = OutputFrame.production
            denudation[:,:,step] = OutputFrame.denudation
            redistribution[:,:,step] = OutputFrame.redistribution
        end
    end
end

end # CaProd
# ~/~ end