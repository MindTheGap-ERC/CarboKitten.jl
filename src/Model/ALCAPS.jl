# ~/~ begin <<docs/src/model-alcap.md#src/Model/ALCAPS.jl>>[init]
module ALCAPS

using Unitful
using ...Config: Box, TimeProperties
using ...SedimentStack: push_sediment!, pop_sediment!

# ~/~ begin <<docs/src/model-alcap.md#alcaps>>[init]
const Rate = typeof(1.0u"m/Myr")
const Amount = typeof(1.0u"m")

@kwdef struct Facies
    viability_range::Tuple{Int,Int}
    activation_range::Tuple{Int,Int}

    maximum_growth_rate::Rate
    extinction_coefficient::typeof(1.0u"m^-1")
    saturation_intensity::typeof(1.0u"W/m^2")

    # though units are in m, this in not an amount.
    # TODO: figure out what this unit means
    # values should be picked rather large, say 10km.
    diffusion_coefficient::typeof(1.0u"m")
end

@kwdef struct Input
    box::Box
    time::TimeProperties

    bedrock_elevation   # (m, m) -> m
    sea_level           # Myr -> m
    subsidence_rate::Rate

    facies::Vector{Facies}
    disintegration_rate::Rate
    insolation::typeof(1.0u"W/m^2")

    sediment_buffer_size::Int
    sediment_buffer_grain::Amount
end
# ~/~ end
# ~/~ begin <<docs/src/model-alcap.md#alcaps>>[1]
struct DisintegrationFrame
    disintegration::Array{Amount,3}    # facies, x, y
end

struct ProductionFrame
    production::Array{Amount,3}        # facies, x, y
end

## Emitting a ModelFrame every iteration allows for
## inspecting the output of an entire run
## For 100x50x3 x 3 x 10000 x 8B ≈ 5GB
struct ModelFrame
    disintegration::Array{Amount,3}    # facies, x, y
    production::Array{Amount,3}
    deposition::Array{Amount,3}
end

mutable struct State
    time::typeof(1.0u"Myr")

    ca::Array{Int}
    ca_priority::Vector{Int}

    sediment_height::Array{Amount,2}   # x, y
    # sediment_buffer stores fractions, so no units
    sediment_buffer::Array{Float64,4}  # facies, x, y, z
end

function initial_state(input)
    height = zeros(Float64, input.box.grid_size...) * u"m"
    for i in CartesianIndices(height)
        height[i] = input.initial_depth(i[1] * input.box.phys_scale)
    end
    n_facies = length(input.facies)
    buffer = zeros(Float64, n_facies, input.box.grid_size..., input.sediment_buffer_size)

    ca = rand(0:n_facies, input.box.grid_size...)
    state = State(0.0u"Myr", ca, 1:n_facies, height)

    step = step_ca(input.box, input.facies)
    for _ = 1:20
        step(state)
    end

    return state
end

function run_model(input)
    state = initial_state(input)
    disintegrate! = disintegration(input)
    produce = production(input)
    transport = transportation(input)

    Channel{ModelFrame}() do ch
        for _ in 1:input.time.steps
            p = produce(state)
            d = disintegrate!(state)

            active_layer = p.production .+ d.disintegration
            sediment = transport(active_layer)

            put!(ch, ModelFrame(d, p, sediment))

            push_sediment!(state.sediment_buffer, sediment ./ input.sediment_buffer_grain |> NoUnits)
            state.sediment_height .+= sum(sediment; dims=1)
            state.time += input.time.𝚫t
        end
    end
end
# ~/~ end

end
# ~/~ end
