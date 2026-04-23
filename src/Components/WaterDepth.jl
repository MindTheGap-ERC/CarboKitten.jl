# ~/~ begin <<docs/src/components/waterdepth.md#src/Components/WaterDepth.jl>>[init]
@compose module WaterDepth
@mixin TimeIntegration, Boxes
using ..Common
using HDF5
using ..TimeIntegration: time, time_axis
using Unitful: ustrip
using CarboKitten.WaveField



export water_depth, SubsidenceModifier

struct SubsidenceModifier
    multiplier::Function   # (x, y, t) -> Float64
end



@kwdef struct Input <: AbstractInput
    sea_level = t -> 0.0u"m"
    initial_topography = (x, y) -> 0.0u"m"
    subsidence_rate = 0.0u"m/Myr"
    subsidence_modifiers::Vector = []
    wave_model::Union{Nothing, WaveField.WaveModel} = nothing
end

function subsidence_rate_matrix(input::AbstractInput)
    sr = input.subsidence_rate
    if sr isa AbstractMatrix
        @assert size(sr) == input.box.grid_size
        return sr
    else
        return fill(sr, input.box.grid_size...)
    end
end

@kwdef mutable struct State <: AbstractState
    sediment_height::Matrix{Height}
    cumulative_subsidence::Matrix{Height} = zeros(Height, 0, 0)
    cumulative_subsidence_hist::Vector{Matrix{Height}} = Matrix{Height}[]
    wdepth_hist::Vector{Matrix{Float64}} = Matrix{Float64}[]
    block_wdepth::Array{Float32,3} = zeros(Float32, 0, 0, 1)
    energy_hist::Vector{Matrix{Float32}} = Matrix{Float32}[]
    block_energy::Array{Float32,3} = zeros(Float32, 0, 0, 1)
end



function initial_state(input::AbstractInput)
    nx, ny = input.box.grid_size

    return State(
        step = 0,
        sediment_height = zeros(Height, nx, ny),
        cumulative_subsidence = zeros(Height, nx, ny),
        cumulative_subsidence_hist = Matrix{Height}[],
        wdepth_hist = Matrix{Float64}[],
        block_wdepth = zeros(Float32, nx, ny, 1),
        energy_hist = Matrix{Float32}[],
        block_energy = zeros(Float32, nx, ny, 1),
    )
end


function initial_topography(input::AbstractInput)
    if input.initial_topography isa AbstractMatrix
        @assert size(input.initial_topography) == input.box.grid_size
        return input.initial_topography
    end

    x, y = box_axes(input.box)
    return input.initial_topography.(x, y')
end


function water_depth(input::AbstractInput)
    eta0 = initial_topography(input)
    sea = input.sea_level
    get_time = time(input)

    return function (state::AbstractState)
        return sea(get_time(state)) .- eta0 .+
               state.cumulative_subsidence .-
               state.sediment_height
    end
end

function write_header(input::AbstractInput, output::AbstractOutput)
    x, y = box_axes(input.box)
    t = time_axis(input)
    set_attribute(output, "initial_topography", initial_topography(input) |> in_units_of(u"m"))
    set_attribute(output, "sea_level", input.sea_level.(t) .|> in_units_of(u"m"))
    set_attribute(output,
              "subsidence_rate",
              ustrip.(u"m/Myr", subsidence_rate_matrix(input)))
end

end
# ~/~ end
