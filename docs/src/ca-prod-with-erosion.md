# Combining CA with production
This model combines BS92 production with the B13 cellular automaton.

## Complete example
This example is running for 10000 steps to 1Myr on a 100 $\times$ 50 grid, starting with a sloped height down to 50m. The `sea_level`, and `initial_depth` arguments are functions. The `phys_scale` argument translate pixels on the grid into physical metres. The `write_interval` indicates to write output every 10 iterations, summing the production over that range. You may copy paste the following code into your own script or notebook, and play around with input values.

``` {.julia .build file=examples/erosion/caps-osc.jl target=data/capes-osc.h5}
using CarboKitten
using CarboKitten.CaProdErosion
using CarboKitten.Burgess2013

#change sinouid function
DEFAULT_INPUT = CarboKitten.CaProdErosion.Input(
    sea_level = t -> 20 * sin(2π * t) + 5 * sin(2π * t / 0.112), 
    subsidence_rate = 50.0,
    initial_depth = x -> x / 2,
    grid_size = (50, 100),
    phys_scale = 1.0,
    Δt = 0.001,
    write_interval = 1,
    time_steps = 2000,
    facies = [
        Burgess2013.Facies((4, 10), (6, 10), 500.0, 0.8, 300, 1000, 2730, 0.5),
        Burgess2013.Facies((4, 10), (6, 10), 400.0, 0.1, 300, 1000, 2730, 0.5),
        Burgess2013.Facies((4, 10), (6, 10), 100.0, 0.005, 300, 1000, 2730, 0.5)
    ],
    insolation = 2000.0,
    temp = 288.0,
    precip = 1000.0,
    pco2 = 10^(-1.5),
    alpha = 2e-3,
    erosion_type = 1

)

CarboKitten.CaProdErosion.main(DEFAULT_INPUT, "data/capes-osc.h5")
```

This writes output to an HDF5 file that you may use for further analysis and visualization.

``` {.julia .build file=examples/erosion/plot-caps-osc.jl target=docs/src/fig/capesosc-crosssection.png deps=data/capes-osc.h5}
module Script
    using CarboKitten.Visualization
    using GLMakie

    function main()
        f = Figure()
        plot_crosssection(f[1,1], "data/capes-osc.h5")
	    save("docs/src/fig/capesosc-crosssection.png", f)
    end
end

Script.main()
```

## Input

- initial depth (function of space)
- sea-level curve (function of time)
- subsidence (function of time)
- facies types

These should all behave as a functions, but could also be some interpolated data. The signs of these quantities should be such that the following equation holds:

$$T + E = S + W$$

Saying Tectonic subsidence plus Eustatic sea-level change equals Sedimentation plus change in Water depth.

``` {.julia #cape-input}
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
```

In the case `write_interval` is not one, we will sum production rates over several iterations of the model before writing to output. In that case sediment production per written frame is no longer limited to a single facies.

## Output

Each iteration of the model, we produce a `Frame`.

``` {.julia #cape-frame}

abstract type Frame end
struct ProductionFrame <: Frame
    production::Array{Float64,3}
end

struct OutputFrame
    production::Array{Float64,3}
    denudation::Array{Float64,2}
    redistribution::Array{Float64,2}
end
```

The frame is used to update a *state* $S$. The frame should be considered a delta for the state. So, we can reproduce the height at each step from the frames.

``` {.julia #cape-state}
# FIXME: deduplicate
mutable struct State
    time::Float64
    height::Array{Float64,2}
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

``` {.julia #cape-model}
# FIXME: deduplicate
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

``` {.julia #cape-model}
# FIXME: deduplicate
function propagator(input::Input)
    <<cape-init-propagator>>
    slopefn = stencil(Float64, Periodic{2}, (3, 3), slope_kernel)
    function (s::State)  # -> Frame
        <<cape-propagate>>
    end
end
```

The propagator keeps the cellular automaton as an internal state, but this may also be considered to be an input function. This may change when you'd want to influence the CA with environmental factors. Then the CA becomes an integral component of the dynamical model. The CA would then have to keep state in the `State` variable. We burn the first 20 iterations of the CA to start with a realistic pattern.

``` {.julia #cape-init-propagator}
n_facies = length(input.facies)
ca_init = rand(0:n_facies, input.grid_size...)
ca = drop(run_ca(Periodic{2}, input.facies, ca_init, 3), 20)

function water_depth(s::State)
    s.height .- input.sea_level(s.time)
end
# prepare functions for erosion
```

Now, to generate a production from a given state, we advance the CA by one step and compute the production accordingly.

``` {.julia #cape-propagate}
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
```

## Updater
Every iteration we update the height variable with the subsidence rate, and add sediments to the height.

``` {.julia #cape-model}
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
```

## Loop

``` {.julia #cape-model}
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
```

``` {.julia file=src/CaProdErosion.jl}
module CaProdErosion

using CarboKitten
using ..Stencil: Periodic, stencil
using ..Utility
#using ..BS92: sealevel_curve
using ..Denudation: DenudationType, denudation, DenudationFrame
using ..Burgess2013

using HDF5
using .Iterators: drop, peel, partition, map, take

<<cape-input>>
<<cape-frame>>
<<cape-state>>
<<cape-model>>

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
```

# Visualizing output

``` {.julia file=src/VisualizationErosion.jl}
module VisualizationErosion

export plot_crosssection

using HDF5
using CairoMakie
using GeometryBasics

function plot_crosssection(pos, datafile)
    # x: 1-d array with x-coordinates
    # t: 1-d array with time-coordinates (n_steps + 1)
    # h[x, t]: height fn, monotonic increasing in time
    # p[x, facies, t]: production rate
    # taken at y = y_max / 2, h[x, 1] is initial height
    n_facies, x, t, h, p = h5open(datafile,"r") do fid
        attr = HDF5.attributes(fid["input"])
        Δt = attr["delta_t"][]
        subsidence_rate = attr["subsidence_rate"][]
        t_end = fid["input/t"][end-1]
        total_subsidence = subsidence_rate * t_end
        total_sediment = sum(fid["sediment"][]; dims=3)
        initial_height = fid["input/height"][]
        center = div(size(total_sediment)[1], 2)
        elevation = cumsum(total_sediment; dims=4)[center,:,1,:] .* Δt .- initial_height .- total_subsidence
        t = fid["input/t"][]
        n_facies = size(fid["sediment"])[3]

        return n_facies,
               fid["input/x"][],
               [t; Δt*attr["time_steps"][]],
               hcat(.- initial_height .- total_subsidence, elevation),
               fid["sediment"][center,:,:,:]
    end

	pts = vec(Point{2,Float64}.(x, h[:,2:end]))
	c = vec(argmax(p; dims=2)[:,1,:] .|> (c -> c[2]))
	rect = Rect2(0.0, 0.0, 1.0, 1.0)
	m_tmp = GeometryBasics.mesh(Tesselation(rect, (100, 1000)))
	m = GeometryBasics.Mesh(pts, faces(m_tmp))

	# pts = vec(Point{2,Float64}.(x, h))
	# c = argmax(p; dims=2)[:,1,:] .|> (c -> c[2])
    # w = size(x)[1]

    # face(idx) = let k = idx[1] + idx[2]*w
    #     TriangleFace(k, k+1, k+1+w), TriangleFace(k+1+w, k+w, k)
    # end

	ax = Axis(pos, xlabel="location", ylabel="depth", limits=((-12,x[end]), nothing))
    # for f in 1:n_facies
    #     locs = CartesianIndices((size(x)[1], size(t)[1] - 1))[c .== f]
    #     triangles = collect(Iterators.flatten(face.(locs)))
    #     m = GeometryBasics.Mesh(pts, triangles)
    #     mesh!(ax, m)
    # end

	mesh!(ax, m, color=c, alpha=0.7)
	for idx in [1,501,1001]
		lines!(ax, x, h[:, idx], color=:black)
		text!(ax, -2.0, h[1, idx]; text="$(t[idx]) Myr", align=(:right, :center))
	end
	for idx in [250,750]
		lines!(ax, x, h[:, idx], color=:black, linewidth=0.5)
	end
end

end
```
