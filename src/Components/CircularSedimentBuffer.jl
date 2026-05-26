# ~/~ begin <<docs/src/components/sediment_buffer.md#src/Components/CircularSedimentBuffer.jl>>[init]
@compose module CircularSedimentBuffer
@mixin Boxes

using ..Common
using CarboKitten.Algorithms.CircularBuffer: push_sediment!, pop_sediment!

export pop_sediment!, push_sediment!

@kwdef struct Input <: AbstractInput
    sediment_buffer_size::Int = 50
    depositional_resolution::Amount = 0.5u"m"
end

@kwdef mutable struct State <: AbstractState
    sediment_buffer::Array{Float64,4}
end

end
# ~/~ end
