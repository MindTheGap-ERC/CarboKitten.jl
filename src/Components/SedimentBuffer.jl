# ~/~ begin <<docs/src/components/sediment_buffer.md#src/Components/SedimentBuffer.jl>>[init]
@compose module SedimentBuffer
@mixin Boxes

using ..Common
using ...Algorithms.CircularBuffers: CircularBufferHost
using ...Interfaces.SedimentBuffers: pop_sediment!, push_sediment!, peek_sediment

@kwdef struct Input <: AbstractInput
    sediment_buffer_size::Int = 50
    depositional_resolution::Amount = 0.5u"m"
end

@kwdef mutable struct State{D, N} <: AbstractState
    sediment_buffer::CircularBufferHost{D, N}
end

end
# ~/~ end
