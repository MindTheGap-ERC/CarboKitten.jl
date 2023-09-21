# ~/~ begin <<docs/src/ca-with-production.md#examples/ca-with-prod.jl>>[init]
module CaProd

using CarboKitten
using CarboKitten.Stencil: Periodic
using CarboKitten.Utility
using CarboKitten.BS92: sealevel_curve
using CarboKitten.Stencil
using CarboKitten.Burgess2013

using .Iterators: drop, peel

# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-input>>[init]
@kwdef struct Input
    sea_level
    subsidence_rate
    initial_depth

    grid_size::NTuple{2,Int}
    phys_scale::Float64
    Δt::Float64

    time_steps::Int
    write_interval::Int

    facies::Vector{Facies}
    insolation::Float64
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-frame>>[init]
struct Frame
    production :: Array{Float64, 3}
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-state>>[init]
mutable struct State
    time :: Float64
    height :: Array{Float64, 2}
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-model>>[init]
function initial_state(input::Input)  # -> State
    height = zeros(Float64, input.grid_size...)
    for i in CartesianIndices(height)
        height[i] = input.initial_depth(i[2] * input.phys_scale)
    end
    return State(0.0, height)
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-model>>[1]
function propagator(input::Input)
    # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-init-propagator>>[init]
    n_facies = length(input.facies)
    ca_init = rand(0:n_facies, input.grid_size...)
    ca = drop(run_ca(Periodic{2}, input.facies, ca_init, 3), 20)

    function water_depth(s::State)
        s.height .- input.sea_level(s.time)
    end
    # ~/~ end
    function (s::State)  # -> Frame
        # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-propagate>>[init]
        result = zeros(Float64, input.grid_size..., n_facies)
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
        s.height .-= sum(Δ.production; dims=3) .* input.Δt
        s.height .+= input.subsidence_rate * input.Δt
        s.time += input.Δt
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

end # CaProd

module Script
    using ..CaProd: Frame, run_model, Input, Facies
    using HDF5
    using .Iterators: partition, map, take
    using CarboKitten.BS92: sealevel_curve

    function stack_frames(fs::Vector{Frame})  # -> Frame
        Frame(sum(f.production for f in fs))
    end

    function main(input::Input)
        x_axis = (0:(input.grid_size[2] - 1)) .* input.phys_scale
        y_axis = (0:(input.grid_size[1] - 1)) .* input.phys_scale
        initial_height = input.initial_depth.(x_axis)
        n_writes = input.time_steps ÷ input.write_interval

        h5open("data/ca-prod.h5", "w") do fid
            gid = create_group(fid, "input")
            gid["x"] = collect(x_axis)
            gid["y"] = collect(y_axis)
            gid["height"] = collect(initial_height)
            gid["t"] = collect((0:(n_writes-1)) .* (input.Δt * input.write_interval))
            attr = attributes(gid)
            attr["delta_t"] = input.Δt
            attr["subsidence_rate"] = input.subsidence_rate

            n_facies = length(input.facies)
            ds = create_dataset(fid, "sediment", datatype(Float64),
                dataspace(input.grid_size..., n_facies, input.time_steps),
                chunk=(input.grid_size..., n_facies, 1))

            results = map(stack_frames, partition(run_model(input), input.write_interval))
            for (step, frame) in enumerate(take(results, n_writes))
                ds[:,:,:,step] = frame.production
            end
        end
    end

    DEFAULT_INPUT = Input(
        sea_level = _ -> 0.0, 
        subsidence_rate = 50.0,
        initial_depth = _ -> 2.0,
        grid_size = (50, 50),
        phys_scale = 1.0,
        Δt = 0.001,
        write_interval = 1,
        time_steps = 1000,
        facies = [
            Facies((4, 10), (6, 10), 500.0, 0.8, 300),
            Facies((4, 10), (6, 10), 400.0, 0.1, 300),
            Facies((4, 10), (6, 10), 100.0, 0.005, 300)
        ],
        insolation = 2000.0
    )
end

Script.main(Script.DEFAULT_INPUT)
# ~/~ end