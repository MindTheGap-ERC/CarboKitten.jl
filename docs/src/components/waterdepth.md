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

export water_depth

@kwdef struct Input <: AbstractInput
    sea_level          # function (t::Time) -> Length
    initial_topography  # function (x::Location, y::Location) -> Length
    subsidence_rate::Rate
end

@kwdef mutable struct State <: AbstractState
    sediment_height::Matrix{Height}
end

function initial_state(input::AbstractInput)
    return State(step=0, sediment_height=zeros(Height, input.box.grid_size...))
end

function water_depth(input::AbstractInput)
    x, y = box_axes(input.box)
    eta0 = input.initial_topography.(x, y')

    return function (state::AbstractState)
        t = TimeIntegration.time(input, state)
        return input.sea_level(t) .- eta0 .+
               (input.subsidence_rate * t) .- state.sediment_height
    end
end

function write_header(fid, input::AbstractInput)
    gid = fid["input"]
    attr = attributes(gid)
    x, y = box_axes(input.box)
    t = TimeIntegration.write_times(input)

    gid["initial_topography"] = input.initial_topography.(x, y') |> in_units_of(u"m")
    gid["sea_level"] = input.sea_level.(t) .|> in_units_of(u"m")
    attr["subsidence_rate"] = input.subsidence_rate |> in_units_of(u"m/Myr")
end

function create_dataset(fid, input::AbstractInput)
    return HDF5.create_dataset(fid, "sediment_height", datatype(Float64),
        dataspace(input.box.grid_size..., input.time.steps),
        chunk=(input.box.grid_size..., 1))
end

end
```
