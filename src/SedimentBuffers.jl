# ~/~ begin <<docs/src/components/sediment_buffer.md#src/SedimentBuffers.jl>>[init]
module SedimentBuffers

export AbstractSedimentBuffer, push_sediment!, pop_sediment!, peek_sediment

using ModuleMixins.Procedural: component
using ..CarboKitten: AbstractState

"""
    AbstractSedimentBuffer

Role for components that store the facies composition of deposited sediment.
An implementing component contributes a `sediment_buffer` field to the `State`.
"""
abstract type AbstractSedimentBuffer end

"""
    push_sediment!(state, sediment::AbstractArray{F, 3})
    push_sediment!(state, i::CartesianIndex{2}, parcel::AbstractVector{F})

Push `sediment` (facies, x, y) onto the sediment buffer, or a `parcel` (facies) onto the column at `i`.
"""
function push_sediment! end

"""
    pop_sediment!(state, amount::AbstractArray{F, 2}, out::AbstractArray{F, 3})
    pop_sediment!(state, i::CartesianIndex{2}, amount::F, parcel::AbstractVector{F})

Pop `amount` (x, y) off the sediment buffer, writing the popped sediment to `out` (facies, x, y). The second form pops a single column, writing to `parcel`.
"""
function pop_sediment! end

"""
    peek_sediment(state, amount::F) -> Array{F, 3}
    peek_sediment(state, i::CartesianIndex{2}, amount::F) -> Vector{F}

Return the sediment that would be popped by `pop_sediment!`, without changing the buffer.
"""
function peek_sediment end

push_sediment!(state::AbstractState, args...) =
    push_sediment!(component(state, AbstractSedimentBuffer), args...)
pop_sediment!(state::AbstractState, args...) =
    pop_sediment!(component(state, AbstractSedimentBuffer), args...)
peek_sediment(state::AbstractState, args...) =
    peek_sediment(component(state, AbstractSedimentBuffer), args...)

end
# ~/~ end
