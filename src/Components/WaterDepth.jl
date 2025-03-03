# ~/~ begin <<docs/src/components/waterdepth.md#src/Components/WaterDepth.jl>>[init]
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
# ~/~ end
