# Water Depth

The `WaterDepth` module computes the water depth, given the bedrock elevation, sea level curve, subsidence rate and current sediment height.

``` {.julia file=src/Components/WaterDepth.jl}
@compose module WaterDepth
@mixin TimeIntegration, Boxes
using ..Common
using HDF5
using CarboKitten.Config: axes

export water_depth

@kwdef struct Input <: AbstractInput
    sea_level          # function (t::Time) -> Length
    bedrock_elevation  # function (x::Location, y::Location) -> Length
    subsidence_rate::Rate
end

@kwdef mutable struct State <: AbstractState
    sediment_height::Matrix{Height}
end

@kwdef struct Frame <: AbstractFrame
    sediment_height::Matrix{Height}
end

function water_depth(input::AbstractInput)
    x, y = axes(input.box)
    eta0 = input.bedrock_elevation.(x, y')

    return function (state::AbstractState)
        t = TimeIntegration.time(input, state)
        return input.sea_level(t) .- eta0 .+
               (input.subsidence_rate * t) .- state.sediment_height
    end
end

function write_header(fid, input::AbstractInput)
    gid = fid["input"]
    attr = attributes(gid)
    x, y = Common.axes(input.box)
    t = (0:n_writes(input)) .* (input.time.Δt * input.time.write_interval)

    gid["bedrock_elevation"] = input.bedrock_elevation.(x, y') |> in_units_of(u"m")
    gid["sea_level"] = input.sea_level.(t) .|> in_units_of(u"m")
    attr["subsidence_rate"] = input.subsidence_rate |> in_units_of(u"m/Myr")
end

function create_dataset(fid, input::AbstractInput)
    return HDF5.create_dataset(fid, "sediment_height", datatype(Float64),
        dataspace(input.box.grid_size..., input.time.steps),
        chunk=(input.box.grid_size..., 1))
end

function write_frame(fid, state::AbstractState, frame::AbstractFrame)
    dataset = fid["sediment_height"]
    dataset[:, :, state.step] = frame.sediment_height
end

end
```

