# Combining CA with production
This model combines BS92 production with the B13 cellular automaton.

## Input

- initial depth (function of space)
- sea-level curve (function of time)
- subsidence (function of time)
- facies types

These should all behave as a functions, but could also be some interpolated data. The signs of these quantities should be such that the following equation holds:

$$T + E = S + W$$

Saying Tectonic subsidence plus Eustatic sea-level change equals Sedimentation plus change in Water depth.

``` {.julia #ca-prod-input}
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
```

In the case `write_interval` is not one, we will sum production rates over several iterations of the model before writing to output. In that case sediment production per written frame is no longer limited to a single facies.

## Output

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

The output is principally all frames produced in the simulation, in a 4-dimensional array. The first two dimensions are x, y positions on the grid, the third is the facies and the fourth dimension is time. We store the output in HDF5, having an `input` group where we store the input data, and a `sediment` dataset containing the aforementioned 4-dimensional output data. Note that these are *production rates*, so to reconstruct the sea floor depth at any time, you need to multiply by $\Delta t * n_w$, where $n_w$ is the `write_interval` and take a cumulative sum.

## Logic

From a dynamical modeling point of view, CarboCAT operates analogous to a forward Euler integration scheme, where some components are actually a discrete model. This means we have one function that generates a `Frame` from a `State`, called the *propagator* $P$ (this is our own nomenclature),

$$P_i: S \to \Delta.$$

The suffix $i$ here is used to indicate that the propagator depends on the input. We'll have a second function $U$ that *updates* the state with the given frame,

$$U: (S, \Delta) \to S.$$

In practice however, the update function changes the state in-place.

## Init

We fill the height map with the initial depth function. It is assumed that the height only depends on the second index.

``` {.julia #ca-prod-model}
function initial_state(input::Input)  # -> State
    height = zeros(Float64, input.grid_size...)
    for i in CartesianIndices(height)
        height[i] = input.initial_depth(i[2] * input.phys_scale)
    end
    return State(0.0, height)
end
```

## Propagator
The propagator computes the production rates (and also erosion) given the state of the model.

``` {.julia #ca-prod-model}
function propagator(input::Input)
    <<ca-prod-init-propagator>>
    function (s::State)  # -> Frame
        <<ca-prod-propagate>>
    end
end
```

The propagator keeps the cellular automaton as an internal state, but this may also be considered to be an input function. This may change when you'd want to influence the CA with environmental factors. Then the CA becomes an integral component of the dynamical model. The CA would then have to keep state in the `State` variable. We burn the first 20 iterations of the CA to start with a realistic pattern.

``` {.julia #ca-prod-init-propagator}
n_facies = length(input.facies)
ca_init = rand(0:n_facies, input.grid_size...)
ca = drop(run_ca(Periodic{2}, input.facies, ca_init, 3), 20)

function water_depth(s::State)
    s.height .- input.sea_level(s.time)
end
```

Now, to generate a production from a given state, we advance the CA by one step and compute the production accordingly.

``` {.julia #ca-prod-propagate}
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
```

## Updater
Every iteration we update the height variable with the subsidence rate, and add sediments to the height.

``` {.julia #ca-prod-model}
function updater(input::Input)
    n_facies = length(input.facies)
    function (s::State, Δ::Frame)
        s.height .-= sum(Δ.production; dims=3) .* input.Δt
        s.height .+= input.subsidence_rate * input.Δt
        s.time += input.Δt
    end
end
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

``` {.julia file=src/CaProd.jl}
module CaProd

export Input, Facies, main

using CarboKitten
using CarboKitten.Stencil: Periodic
using CarboKitten.Utility
using CarboKitten.BS92: sealevel_curve
using CarboKitten.Stencil
using CarboKitten.Burgess2013

using HDF5
using .Iterators: drop, peel, partition, map, take

<<ca-prod-input>>
<<ca-prod-frame>>
<<ca-prod-state>>
<<ca-prod-model>>

function stack_frames(fs::Vector{Frame})  # -> Frame
    Frame(sum(f.production for f in fs))
end

function main(input::Input, output::String)
    x_axis = (0:(input.grid_size[2] - 1)) .* input.phys_scale
    y_axis = (0:(input.grid_size[1] - 1)) .* input.phys_scale
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

end # CaProd
```

## Case 1
The first case uses the same settings as Burgess 2013: an initial depth of 2m, subsidence rate of 50 m/Myr and constant sea level.

``` {.julia file=examples/ca-with-prod.jl}
using CarboKitten.CaProd

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

main(DEFAULT_INPUT, "data/ca-prod.h5")
```

## Case 2
For the second case, we start with a slope.

``` {.julia file=examples/cap-slope.jl}
using CarboKitten.CaProd

DEFAULT_INPUT = Input(
    sea_level = _ -> 0.0, 
    subsidence_rate = 50.0,
    initial_depth = x -> x,
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

main(DEFAULT_INPUT, "data/ca-prod-slope.h5")
```

# Visualizing output

![](crosssection.png)

```@example
using HDF5
using CairoMakie
using GeometryBasics

function main()
    x, h, p = h5open("../../data/ca-prod-slope.h5","r") do fid
        attr = HDF5.attributes(fid["input"])
        Δt = attr["delta_t"][]
        subsidence_rate = attr["subsidence_rate"][]
        t_end = fid["input/t"][end-1]
        total_subsidence = subsidence_rate * t_end
        total_sediment = sum(fid["sediment"][]; dims=3)
        initial_height = fid["input/height"][]
        # compute cumulative production and correct for initial height and subsidence
        elevation = cumsum(total_sediment; dims=4)[1,:,1,:] .* Δt .- initial_height .- total_subsidence
        return fid["input/x"][], elevation, fid["sediment"][1,:,:,:]
    end

    pts = vec(Point{2,Float64}.(x, h))
    # colours determined by which facies has maximum production
    c = vec(argmax(p; dims=2)[:,1,:] .|> (c -> c[2]))
    # create mesh for regular grid
    rect = Rect2(0.0, 0.0, 1.0, 1.0)
    m_tmp = GeometryBasics.mesh(Tesselation(rect, (50, 1000)))
    # and move vertices to x, h positions
    m = GeometryBasics.Mesh(pts, faces(m_tmp))

    f = Figure()
    ax = Axis(f[1, 1], xlabel="location", ylabel="depth")
    mesh!(ax, m, color=c, alpha=0.7)
    save("crosssection.png", f)
end

main()
```
