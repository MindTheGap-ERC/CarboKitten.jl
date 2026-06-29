# Water Depth

```component-dag
CarboKitten.Components.WaterDepth
```

The `WaterDepth` module computes the water depth, given the bedrock elevation, sea level curve, subsidence rate and current sediment height.

## Input

- `initial_topography(x, y)` (a.k.a. initial depth) should be a function taking two coordinates in units of meters, returning an elevation also in meters.
- `sea_level(t)` should be a function taking a time in millions of years (Myr) returning the eustatic sealevel. This could also be an interpolated table.
- `subsidence_rate` a constant rate of subsidence in m/Myr.

The signs of these quantities should be such that the following equation holds:

$$T + E = S + W,$$

saying Tectonic subsidence plus Eustatic sea-level change equals Sedimentation plus change in Water depth.

``` {.julia file=src/Components/WaterDepth.jl}
@compose module WaterDepth
@mixin TimeIntegration, Boxes
using ..Common
using HDF5
using ..TimeIntegration: time, time_axis

export water_depth, subsider, initial_topography

@kwdef struct Input <: AbstractInput
    sea_level = t -> 0.0u"m"
    initial_topography = (x, y) -> 0.0u"m"
    subsidence_rate::Rate = 0.0u"m/Myr"
end

@kwdef mutable struct State <: AbstractState
    bathymetry::Matrix{Height}
end

function initial_state(input::AbstractInput)
    bathymetry = initial_topography(input)
    return State(step=0, bathymetry=bathymetry)
end

function initial_topography(input::AbstractInput)
    if input.initial_topography isa AbstractMatrix
        @assert size(input.initial_topography) == input.box.grid_size
        return input.initial_topography
    end

    x, y = box_axes(input.box)
    return input.initial_topography.(x, y')
end

function subsider(input::AbstractInput)
    Δσ = input.subsidence_rate * input.time.Δt

    function (state::AbstractState)
        state.bathymetry .-= Δσ
    end
end

function water_depth(input::AbstractInput)
    sea_level = input.sea_level
    get_time = time(input)

    return function (state::AbstractState)
        t = get_time(state)
        return sea_level(t) .- state.bathymetry
    end
end

function write_header(input::AbstractInput, output::AbstractOutput)
    x, y = box_axes(input.box)
    t = time_axis(input)
    set_attribute(output, "initial_topography", initial_topography(input) |> in_units_of(u"m"))
    set_attribute(output, "sea_level", input.sea_level.(t) .|> in_units_of(u"m"))
    set_attribute(output, "subsidence_rate", input.subsidence_rate |> in_units_of(u"m/Myr"))
end

end
```
