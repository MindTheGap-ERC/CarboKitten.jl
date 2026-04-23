# ~/~ begin <<docs/src/components/sediment_buffer.md#src/Components/SedimentBuffer.jl>>[init]
@compose module SedimentBuffer
@mixin Boxes

using ..Common
using CarboKitten.SedimentStack: pop_sediment!, push_sediment!, peek_sediment
using Random

export pop_sediment!, push_sediment!, peek_sediment

@kwdef struct Input <: AbstractInput
    sediment_buffer_size::Int = 50
    depositional_resolution::Amount = 0.5u"m"
end

@kwdef mutable struct State <: AbstractState
    sediment_buffer::Array{Float64,4}

    # --- in-memory facies block accumulators ---
    block_cube::Array{UInt8,3} = zeros(UInt8, 0, 0, 1)
    block_topk::Array{Int,2} = zeros(Int, 0, 0)

    # --- Wheeler / time-history storage ---
    production_hist::Vector{Array{Float64,3}} = Array{Float64,3}[]
    disintegration_hist::Vector{Array{Float64,3}} = Array{Float64,3}[]
    deposition_hist::Vector{Array{Float64,3}} = Array{Float64,3}[]
    sediment_thickness_hist::Vector{Matrix{Float64}} = Matrix{Float64}[]
    time_hist::Vector{Float64} = Float64[]
end

end

# ~/~ end