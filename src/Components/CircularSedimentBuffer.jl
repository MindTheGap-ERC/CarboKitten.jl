# ~/~ begin <<docs/src/components/sediment_buffer.md#src/Components/CircularSedimentBuffer.jl>>[init]
@component module CircularSedimentBuffer
@mixin Boxes, FaciesBase
@implements AbstractSedimentBuffer

using ..Common
using ..FaciesBase: n_facies
using ...SedimentBuffers: AbstractSedimentBuffer
import ...SedimentBuffers: push_sediment!, pop_sediment!, peek_sediment
using ...Algorithms: CircularBuffers
using ...Algorithms.CircularBuffers: CircularBuffer

struct Input <: InputPart end

@kwdef mutable struct State <: StatePart
    sediment_buffer::CircularBuffer{Float64}
end

Common.initial_state(::Input, input::AbstractInput) = State(
    sediment_buffer=CircularBuffer{Float64}(input.sediment_buffer_size, n_facies(input), input.box.grid_size))

push_sediment!(s::State, args...) = CircularBuffers.push_sediment!(s.sediment_buffer, args...)
pop_sediment!(s::State, args...) = CircularBuffers.pop_sediment!(s.sediment_buffer, args...)
peek_sediment(s::State, args...) = CircularBuffers.peek_sediment(s.sediment_buffer, args...)

end
# ~/~ end
