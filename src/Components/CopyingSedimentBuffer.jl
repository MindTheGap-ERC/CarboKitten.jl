# ~/~ begin <<docs/src/components/sediment_buffer.md#src/Components/CopyingSedimentBuffer.jl>>[init]
@component module CopyingSedimentBuffer
@mixin Boxes, FaciesBase
@implements AbstractSedimentBuffer

using ..Common
using ..FaciesBase: n_facies
using ...SedimentBuffers: AbstractSedimentBuffer
import ...SedimentBuffers: initial_sediment_buffer, push_sediment!, pop_sediment!, peek_sediment
import ModuleMixins.Procedural: default_implementation
using ...Algorithms: CopyingBuffers

struct Input <: InputPart end

@kwdef mutable struct State <: StatePart
    sediment_buffer::Array{Float64, 4}
end

initial_sediment_buffer(::Input, input::AbstractInput) =
    zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...)

push_sediment!(s::State, args...) = CopyingBuffers.push_sediment!(s.sediment_buffer, args...)
pop_sediment!(s::State, args...) = CopyingBuffers.pop_sediment!(s.sediment_buffer, args...)
peek_sediment(s::State, args...) = CopyingBuffers.peek_sediment(s.sediment_buffer, args...)

default_implementation(::Type{AbstractSedimentBuffer}) = CopyingSedimentBuffer

end
# ~/~ end
