# Combining CA with production
This model combines BS92 production with the B13 cellular automaton.

We need as an input:

- initial depth (function of space)
- sea-level curve (function of time)
- subsidence (function of time)

These should all behave as a functions, but could also be some interpolated data. The signs of these quantities should be such that the following equation holds:

$$T + E = S + W$$

Saying Tectonic subsidence plus Eustatic sea-level change equals Sedimentation plus change in Water depth.

``` {.julia #ca-prod-input}
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
```

Each iteration of the model, we produce a `Frame`.

``` {.julia #ca-prod-frame}
struct Frame
    production :: Array{Float64, 3}
end
```

The frame is used to update a *state* $S$. The frame should be considered a delta for the state. So, we can reproduce the height at each step from the frames.

``` {.julia #ca-prod-state}
mutable struct State
    time :: Float64
    height :: Array{Float64, 2}
end
```

From a dynamical modeling point of view, CarboCAT operates analogous to a forward Euler integration scheme, where some components are actually a discrete model. This means we have one function that generates a `Frame` from a `State`, called the *propagator* $P$ (this is our own nomenclature),

$$P_i: S \to \Delta.$$

The suffix $i$ here is used to indicate that the propagator depends on the input.

``` {.julia #ca-prod-model}
function initial_state(input::Input)  # -> State
    <<ca-prod-init>>
end

function propagator(input::Input)
    <<ca-prod-init-propagator>>
    function (s::State)  # -> Frame
        <<ca-prod-propagate>>
    end
end
```

We'll have a second function $U$ that *updates* the state with th given frame,

$$U: (S, \Delta) \to S.$$

In practice however, the update function changes the state in-place.

``` {.julia #ca-prod-model}
function updater(input::Input)
    n_facies = length(input.facies)
    function (s::State, Δ::Frame)
        <<ca-prod-update>>
    end
end
```

## Init

We fill the height map with the initial depth function. It is assumed that the height only depends on the second index.

``` {.julia #ca-prod-init}
# n_facies = length(input.facies)
height = zeros(Float64, input.grid_size...)
for i in CartesianIndices(height)
    height[i] = input.initial_depth(i[2] * input.phys_scale)
end
return State(0.0, height)
```

## Propagator
The propagator keeps the cellular automaton as an internal state, but this may also be considered to be an input function. This may change when you'd want to influence the CA with environmental factors. Then the CA becomes an integral component of the dynamical model. The CA would then have to keep state in the `State` variable. We burn the first 20 iterations of the CA to start with a realistic pattern.

``` {.julia #ca-prod-init-propagator}
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
```

Now, to generate a production from a given state, we advance the CA by one step and compute the production accordingly.

``` {.julia #ca-prod-propagate}
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
```

## Updater
Every iteration we update the height variable with the subsidence rate, and add sediments to the height.

``` {.julia #ca-prod-update}
s.height .-= sum(Δ.production; dims=3) .* input.Δt
s.height .+= input.subsidence_rate * input.Δt
s.time += input.Δt
```

## Loop

``` {.julia #ca-prod-model}
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
```

``` {.julia file=examples/ca-with-prod.jl}
module CaProd

using CarboKitten
using CarboKitten.Stencil: Periodic
using CarboKitten.Utility
using CarboKitten.BS92: sealevel_curve
using CarboKitten.Stencil

using .Iterators: drop, peel

<<ca-prod-input>>
<<ca-prod-frame>>
<<ca-prod-state>>


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

<<ca-prod-model>>
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
```

# Visualizing output

``` {.julia file=examples/visualize_production.jl}
module Script
    using HDF5
    using Plots

    function main()
        h5open("data/ca-prod.h5", "r") do fid
            total_sediment = sum(fid["sediment"][:,:,:,:]; dims=3)
            initial_height = fid["input"]["height"][:]
            height = initial_height' .- cumsum(total_sediment; dims=4)
            height[1,:,1,:]
        end
    end
end

Script.main()
```
