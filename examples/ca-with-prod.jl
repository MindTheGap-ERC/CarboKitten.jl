# ~/~ begin <<docs/src/ca-with-production.md#examples/ca-with-prod.jl>>[init]
module CaProd

using CarboKitten
using CarboKitten.Stencil: Periodic
using CarboKitten.Utility
using CarboKitten.BS92: sealevel_curve
using CarboKitten.Stencil

using .Iterators: drop, peel

# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-input>>[init]
struct Facies
    viability_range::Tuple{Int, Int}
    activation_range::Tuple{Int, Int}

    maximum_growth_rate::Float64
    extinction_coefficient::Float64
    saturation_intensity::Float64
end

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

function production(insolation::Float64, facies::Facies, water_depth::Float64)
    gₘ = facies.maximum_growth_rate
    I₀ = insolation
    Iₖ = facies.saturation_intensity
    w = water_depth
    k = facies.extinction_coefficient
    return w > 0.0 ? gₘ * tanh(I₀/Iₖ * exp(-w * k)) : 0.0
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


cycle_permutation(n_species::Int) =
    (circshift(1:n_species, x) for x in Iterators.countfrom(0))

function rules(facies::Vector{Facies})
    function (neighbourhood::Matrix{Int}, order::Vector{Int})
        cell_facies = neighbourhood[3, 3]
        neighbour_count(f) = sum(neighbourhood .== f)
        if cell_facies == 0
            for f in order
                n = neighbour_count(f)
                (a, b) = facies[f].activation_range
                if a <= n && n <= b
                    return f
                end
            end
            0
        else
            n = neighbour_count(cell_facies) - 1
            (a, b) = facies[cell_facies].viability_range
            (a <= n && n <= b ? cell_facies : 0)
        end
    end    
end

function run_ca(::Type{B}, facies::Vector{Facies}, init::Matrix{Int}, n_species::Int) where {B <: Boundary{2}}
    r = rules(facies)
    Channel{Matrix{Int}}() do ch
        target = Matrix{Int}(undef, size(init))
        put!(ch, init)
        stencil_op = stencil(Int, B, (5, 5), r)
        for perm in cycle_permutation(n_species)
            stencil_op(init, target, perm)
            init, target = target, init
            put!(ch, init)
        end
    end
end

# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-model>>[init]
function initial_state(input::Input)  # -> State
    # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-init>>[init]
    # n_facies = length(input.facies)
    height = zeros(Float64, input.grid_size...)
    for i in CartesianIndices(height)
        height[i] = input.initial_depth(i[2] * input.phys_scale)
    end
    return State(0.0, height)
    # ~/~ end
end

function propagator(input::Input)
    # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-init-propagator>>[init]
    n_facies = length(input.facies)
    ca_init = rand(0:n_facies, input.grid_size...)
    # ca = drop(run_ca(Periodic{2}, input.facies, ca_init, 3), 20)
    ca = run_ca(Periodic{2}, input.facies, ca_init, 3)
    for _ in 1:20
        take!(ca)
    end

    function water_depth(s::State)
        s.height .- input.sea_level(s.time)
    end
    # ~/~ end
    function (s::State)  # -> Frame
        # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-propagate>>[init]
        result = zeros(Float64, input.grid_size..., n_facies)
        # facies_map, ca = peel(ca)
        facies_map = take!(ca)
        w = water_depth(s)
        for idx in CartesianIndices(facies_map)
            f = facies_map[idx]
            if f == 0
                continue
            end
            result[Tuple(idx)..., f] = production(input.insolation, input.facies[f], w[idx])
        end
        return Frame(result)
        # ~/~ end
    end
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-model>>[1]
function updater(input::Input)
    n_facies = length(input.facies)
    function (s::State, Δ::Frame)
        # ~/~ begin <<docs/src/ca-with-production.md#ca-prod-update>>[init]
        s.height .-= sum(Δ.production; dims=3) .* input.Δt
        s.height .+= input.subsidence_rate * input.Δt
        s.time += input.Δt
        # ~/~ end
    end
end
# ~/~ end
# ~/~ begin <<docs/src/ca-with-production.md#ca-prod-model>>[2]
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
        sea_level = sealevel_curve(),
        subsidence_rate = 0.0,
        initial_depth = x -> x,
        grid_size = (50, 50),
        phys_scale = 3.0,
        Δt = 10,
        write_interval = 100,
        time_steps = 8000,
        facies = [
            Facies((4, 10), (6, 10), 0.05, 0.8, 300),
            Facies((4, 10), (6, 10), 0.04, 0.1, 300),
            Facies((4, 10), (6, 10), 0.01, 0.005, 300)
        ],
        insolation = 2000.0
    )
end

Script.main(Script.DEFAULT_INPUT)
# ~/~ end