# ~/~ begin <<docs/src/ca-with-production.md#src/CaProd.jl>>[init]
module CaProd

using CarboKitten
using ..Utility
using ..Config: Box, TimeProperties
using ..BoundaryTrait: Periodic
using ..Utility
# using CarboKitten.BS92: sealevel_curve
using ..Stencil
using ..Burgess2013

using HDF5
using Unitful
using .Iterators: drop, peel, partition, map, take

# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-input>>[init]
@kwdef struct Input
    box :: Box
    time :: TimeProperties

    sea_level       # Myr -> m
    subsidence_rate::typeof(1.0u"m/Myr")
    initial_depth   # m -> m

    facies::Vector{Facies}
    insolation::typeof(1.0u"W/m^2")
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-frame>>[init]
struct Frame
    production::Array{typeof(1.0u"m/Myr"),3}
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-state>>[init]
mutable struct State
    time::typeof(1.0u"Myr")
    height::Array{typeof(1.0u"m"),2}
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-model>>[init]
function initial_state(input)  # -> State
    height = zeros(Float64, input.box.grid_size...) * u"m"
    for i in CartesianIndices(height)
        height[i] = input.initial_depth(i[1] * input.box.phys_scale)
    end
    return State(0.0u"Myr", height)
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-model>>[1]
function propagator(input)
    # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-init-propagator>>[init]
    n_facies = length(input.facies)
    ca_init = rand(0:n_facies, input.box.grid_size...)
    ca = drop(run_ca(Periodic{2}, input.facies, ca_init, 3), 20)

    function water_depth(s)
        s.height .- input.sea_level(s.time)
    end
    # ~/~ end
    function (s)  # -> Frame
        # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-propagate>>[init]
        result = zeros(typeof(0.0u"m/Myr"), input.box.grid_size..., n_facies)
        facies_map, ca = peel(ca)
        w = water_depth(s)
        for idx in CartesianIndices(facies_map)
            f = facies_map[idx]
            if f == 0
                continue
            end
            result[Tuple(idx)..., f] = production_rate(input.insolation, input.facies[f], w[idx])
        end
        return Frame(result)
        # ~/~ end
    end
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-model>>[2]
function updater(input::Input)
    n_facies = length(input.facies)
    function (s::State, Δ::Frame)
        s.height .-= sum(Δ.production; dims=3) .* input.time.Δt
        s.height .+= input.subsidence_rate * input.time.Δt
        s.time += input.time.Δt
    end
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-model>>[3]
function run_model(input::Input)
    Channel{Frame}() do ch
        s = initial_state(input)
        p = propagator(input)
        u = updater(input)

        while true
            Δ = p(s)
            put!(ch, Δ)
            u(s, Δ)
        end
    end
end
# ~/~ end

function stack_frames(fs::Vector{Frame})  # -> Frame
    Frame(sum(f.production for f in fs))
end

function main(input::Input, output::String)
    y_axis = (0:(input.box.grid_size[2]-1)) .* input.box.phys_scale
    x_axis = (0:(input.box.grid_size[1]-1)) .* input.box.phys_scale
    initial_height = input.initial_depth.(x_axis)
    n_writes = input.time.steps ÷ input.time.write_interval

    h5open(output, "w") do fid
        gid = create_group(fid, "input")
        gid["x"] = collect(x_axis) |> in_units_of(u"m")
        gid["y"] = collect(y_axis) |> in_units_of(u"m")
        gid["height"] = collect(initial_height) |> in_units_of(u"m")
        gid["t"] = collect((0:(n_writes-1)) .* (input.time.Δt * input.time.write_interval)) |> in_units_of(u"Myr")
        attr = attributes(gid)
        attr["delta_t"] = input.time.Δt |> in_units_of(u"Myr")
        attr["write_interval"] = input.time.write_interval
        attr["time_steps"] = input.time.steps
        attr["subsidence_rate"] = input.subsidence_rate |> in_units_of(u"m/Myr")

        n_facies = length(input.facies)
        ds = create_dataset(fid, "sediment", datatype(Float64),
            dataspace(input.box.grid_size..., n_facies, input.time.steps),
            chunk=(input.box.grid_size..., n_facies, 1))

        results = map(stack_frames, partition(run_model(input), input.time.write_interval))
        for (step, frame) in enumerate(take(results, n_writes))
            ds[:, :, :, step] = frame.production |> in_units_of(u"m/Myr")
        end
    end
end

end # CaProd
# ~/~ end