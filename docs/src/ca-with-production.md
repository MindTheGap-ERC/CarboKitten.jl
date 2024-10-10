# Combining CA with production

This model combines BS92 production with the B13 cellular automaton.

## Complete example

This example is running for 10000 steps to 1Myr on a 100 $\times$ 50 grid, starting with a sloped height down to 50m. The `sea_level`, and `initial_depth` arguments are functions. The `phys_scale` argument translate pixels on the grid into physical metres. The `write_interval` indicates to write output every 10 iterations, summing the production over that range. You may copy paste the following code into your own script or notebook, and play around with input values.

``` {.julia .task file=examples/production-only/caps-osc.jl}
#| creates: data/caps-osc.h5
#| requires: src/CaProd.jl

module Script
using CarboKitten.BoundaryTrait: Shelf
using CarboKitten.Config: Box, TimeProperties
using CarboKitten.Burgess2013.Config: MODEL1
using CarboKitten.CaProd
using Unitful

const PERIOD = 200.0u"kyr"
const AMPLITUDE = 4.0u"m"

const DEFAULT_INPUT = CaProd.Input(
  box = Box{Shelf}(
    grid_size = (100, 50),
    phys_scale = 0.15u"km"
    # equivalent:
    # phys_scale = 150"m"
  ),
  time = TimeProperties(
    Δt = 0.0001u"Myr",
    # equivalent: 
    # Δt = 1u"kyr",
    steps = 10000,
    write_interval = 10
  ),
  sea_level = t -> AMPLITUDE * sin(2π * t / PERIOD), 
  subsidence_rate=50.0u"m/Myr",
  initial_depth=x -> x / 300.0,
  facies=MODEL1,
  insolation=400.0u"W/m^2"
)
end

Script.CaProd.main(Script.DEFAULT_INPUT, "data/caps-osc.h5")
```

This writes output to an HDF5 file that you may use for further analysis and visualization.

``` {.julia .task file=examples/plot-caps-osc.jl}
#| creates: docs/src/_fig/b13-capsosc-crosssection.png
#| requires: data/caps-osc.h5 ext/VisualizationExt.jl
#| collect: figures

module Script
    using CairoMakie
    using GeometryBasics
    using CarboKitten.Visualization

    function main()
        f = Figure()
        plot_crosssection(f[1,1], "data/caps-osc.h5")
       save("docs/src/_fig/b13-capsosc-crosssection.png", f)
    end
end

Script.main()
```

![Stratigraphy, production and subsidence under oscillating sea level.](fig/b13-capsosc-crosssection.png)

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
    box :: Box
    time :: TimeProperties

    sea_level       # Myr -> m
    subsidence_rate::typeof(1.0u"m/Myr")
    initial_depth   # m -> m

    facies::Vector{Facies}
    insolation::typeof(1.0u"W/m^2")
end
```

In the case `write_interval` is not one, we will sum production rates over several iterations of the model before writing to output. In that case sediment production per written frame is no longer limited to a single facies.

## Output

Each iteration of the model, we produce a `Frame`.

``` {.julia #ca-prod-frame}
struct Frame
    production::Array{typeof(1.0u"m/Myr"),3}
end
```

The frame is used to update a *state* $S$. The frame should be considered a delta for the state. So, we can reproduce the height at each step from the frames.

``` {.julia #ca-prod-state}
mutable struct State
    time::typeof(1.0u"Myr")
    ca::Array{Int}
    ca_priority::Vector{Int}
    height::Array{typeof(1.0u"m"),2}
end
```

The output is principally all frames produced in the simulation, in a 4-dimensional array. The first two dimensions are x, y positions on the grid, the third is the facies and the fourth dimension is time. We store the output in HDF5, having an `input` group where we store the input data, and a `sediment` dataset containing the aforementioned 4-dimensional output data. Note that these are *production rates*, so to reconstruct the sea floor depth at any time, you need to multiply by $\Delta t * n_w$, where $n_w$ is the `write_interval` and take a cumulative sum.

## Logic

From a dynamical modeling point of view, CarboCAT operates analogous to a forward Euler integration scheme, where some components are actually a discrete model. This means we have one function that generates a `Frame` from a `State`, called the *propagator* $P$ (this is our own nomenclature),

$$P_i: S \to \Delta.$$

The suffix $i$ here is used to indicate that the propagator depends on the input. We'll have a second function $U$ that *updates* the state with the given frame,

$$U: (S, \Delta) \to S.$$

In practice however, the update function changes the state in-place.

## Implementation

``` {.julia file=src/Model/CAP.jl}
@compose module CAP
@mixin Tag, H5Writer, CAProduction

using ..Common
using ..CAProduction: production
using ..TimeIntegration
using ..WaterDepth
using ModuleMixins: @for_each

export Input, Facies, run

function initial_state(input::Input)
    ca_state = CellularAutomaton.initial_state(input)
    for _ in 1:20
        CellularAutomaton.stepper(input)(ca_state)
    end

    sediment_height = zeros(Height, input.box.grid_size...)
    return State(
        step=0, sediment_height=sediment_height,
        ca=ca_state.ca, ca_priority=ca_state.ca_priority)
end

function step!(input::Input)
    τ = production(input)
    step_ca = CellularAutomaton.stepper(input)

    function (state::State)
        if mod(state.step, input.ca_interval) == 0
            step_ca(state)
        end

        prod = τ(state)
        Δη = sum(prod; dims=1)[1, :, :]
        state.sediment_height .+= Δη
        state.step += 1

        return H5Writer.DataFrame(
            production = prod,
            deposition = prod)
    end
end

function write_header(fid, input::AbstractInput)
    @for_each(P -> P.write_header(fid, input), PARENTS)
end
end
```

## Case 1

The first case uses the same settings as Burgess 2013: an initial depth of 2m, subsidence rate of 50 m/Myr and constant sea level.

``` {.julia .task file=examples/production-only/ca-uniform.jl}
#| creates: data/ca-prod.h5
#| requires: src/CaProd.jl

module Script
using CarboKitten.CaProd
using CarboKitten.BoundaryTrait
using CarboKitten.Config: Box, TimeProperties
using CarboKitten.Burgess2013.Config: MODEL1
using Unitful

const DEFAULT_INPUT = CaProd.Input(
  box = Box{Periodic{2}}(
    grid_size = (50, 50),
    phys_scale = 1.0u"m"
  ),
  time = TimeProperties(
    Δt = 0.001u"Myr",
    steps = 1000,
    write_interval = 1
  ),
  sea_level=_ -> 0.0u"m",
  subsidence_rate=50.0u"m/Myr",
  initial_depth=_ -> 2.0u"m",
  facies=MODEL1,
  insolation=2000.0u"W/m^2"
)
end

CaProd.main(Script.DEFAULT_INPUT, "data/ca-prod.h5")
```

## Case 2

For the second case, we start with a slope.

``` {.julia .task file=examples/production-only/ca-slope.jl}
#| creates: data/ca-prod-slope.h5
#| requires: src/CaProd.jl

module Script
using CarboKitten.BoundaryTrait: Shelf
using CarboKitten.Config: Box, TimeProperties
using CarboKitten.Burgess2013.Config: MODEL1
using CarboKitten.CaProd
using Unitful

const DEFAULT_INPUT = CaProd.Input(
  box = Box{Shelf}(
    grid_size = (100, 50),
    phys_scale = 0.15u"km"
  ),
  time = TimeProperties(
    Δt = 0.001u"Myr",
    steps = 1000,
    write_interval = 1
  ),
  sea_level=_ -> 0.0u"m",
  subsidence_rate=50.0u"m/Myr",
  initial_depth=x -> x / 300.0,
  facies=MODEL1,
  insolation=400.0u"W/m^2"
)
end

Script.CaProd.main(Script.DEFAULT_INPUT, "data/ca-prod-slope.h5")
```

``` {.julia .task file=examples/production-only/plot-cap-slope.jl}
#| creates: docs/src/_fig/b13-crosssection.png
#| requires: data/ca-prod-slope.h5 ext/VisualizationExt.jl
#| collect: figures

module Script
using CairoMakie
using CairoMakie.Export: read_slice
using CarboKitten.Visualization

function main()
    header, data = read_slice("data/ca-prod-slope.h5", 25)
    fig = sediment_profile(header, data)
    save("docs/src/_fig/b13-crosssection.png", fig)
end
end

Script.main()
```

![Stratigraphy; production and subsidence.](fig/b13-crosssection.png)
