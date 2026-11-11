# Sediment Buffers

```component-dag
CarboKitten.Components.SedimentBuffer
```

For our models of transport and denudation it is important to remember the facies of the sediment for some time into the past. One way to do this, is to remember all contributions of sediment in a stack. Every time we transport or erode sediment, we can pop parcels from this stack. In a three-dimensional model, we need a 2d grid of stacks. Each stack would have its own memory management (which is computationally expensive), and most resources are spent on areas with very little accretion. (In fact, this is what CarboCAT does. We believe this design choice is the main contributor to the difference in run-time between CarboCAT and CarboKitten).

Instead, we choose a fixed size sediment buffer. Each cell in the buffer represents a parcel of sediment, where we store the relative fractions of each contributing facies. This buffer is only used to determine the facies of disintegrated sediment. The output of the overal model is still the amount of sediment for each iteration.

We can't stress this enough: any inaccuracy in using a fixed size buffer with a chosen granularity only impacts the precision of the composition of transported sediment. Even then, the schema is conservative: no sediment is lost unless erosion is so rampant that it eats through the entire sediment stack. In that case, a simulation should be run with a larger buffer.

There are two implementations of the sediment buffer:

- The [copying stack](../algorithms/copying_stack.md) keeps the top of the stack at the first row, copying rows down or up when the top overflows or is emptied.
- The [circular stack](../algorithms/circular_stack.md) remembers the position of the top of the stack, so that no copying is needed.

## Role

Components that need a sediment buffer don't need to know which implementation is used. The `SedimentBuffer` component requires the `AbstractSedimentBuffer` role, which is implemented by the `CopyingSedimentBuffer` (the default) and `CircularSedimentBuffer` components. Each of these contributes a `sediment_buffer` to the `State`, and implements the following methods for its own `State` part. When called on a complete `State`, these methods forward to the part that implements the role.

All amounts are unit-free, measured in units of the depositional resolution.

``` {.julia file=src/SedimentBuffers.jl}
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
```

## Component

The `SedimentBuffer` component holds the parameters of the buffer and keeps track of the total `sediment_thickness`. The `push_sediment` and `pop_sediment` functions convert sediment amounts to units of depositional resolution and update the bathymetry.

``` {.julia file=src/Components/SedimentBuffer.jl}
@component module SedimentBuffer
@mixin Boxes, FaciesBase, WaterDepth
@requires AbstractSedimentBuffer

using Unitful

using ..Common
using ..FaciesBase: n_facies
using ...SedimentBuffers: AbstractSedimentBuffer, push_sediment!, pop_sediment!

export pop_sediment, push_sediment

@kwdef struct Input <: InputPart
    sediment_buffer_size::Int = 50
    depositional_resolution::Amount = 0.5u"m"
end

@kwdef mutable struct State <: StatePart
    sediment_thickness::Array{Height, 2}
end

Common.initial_state(::Input, input::AbstractInput) =
    State(sediment_thickness=zeros(Height, input.box.grid_size...))

function push_sediment(input::AbstractInput)
    res = input.depositional_resolution
    n_f = n_facies(input)
    n_g = input.box.grid_size
    parcel = Vector{Float64}(undef, n_f)

    function (state::AbstractState, sediment::Array{Amount, 3})
        for i in CartesianIndices(n_g)
            total = sum(@view sediment[:, i[1], i[2]])
            state.sediment_thickness[i] += total
            state.bathymetry[i] += total
            parcel .= @view(sediment[:, i[1], i[2]]) ./ res .|> NoUnits
            push_sediment!(state, i, parcel)
        end
    end
end

function pop_sediment(input::AbstractInput)
    res = input.depositional_resolution
    n_f = n_facies(input)
    n_g = input.box.grid_size
    parcel = Vector{Float64}(undef, n_f)

    function (state::AbstractState, amount::Array{Amount, 2}, out::Array{Amount, 3})
        for i in CartesianIndices(n_g)
            state.sediment_thickness[i] -= amount[i]
            state.bathymetry[i] -= amount[i]
            pop_sediment!(state, i, amount[i] / res |> NoUnits, parcel)
            view(out, :, i[1], i[2]) .= parcel .* res
        end
    end
end

end
```

## Implementations

Both implementations are thin wrappers around their algorithms. The `CopyingSedimentBuffer` is used when no other implementation is selected.

``` {.julia file=src/Components/CopyingSedimentBuffer.jl}
@component module CopyingSedimentBuffer
@mixin Boxes, FaciesBase
@implements AbstractSedimentBuffer

using ..Common
using ..FaciesBase: n_facies
using ...SedimentBuffers: AbstractSedimentBuffer
import ...SedimentBuffers: push_sediment!, pop_sediment!, peek_sediment
import ModuleMixins.Procedural: default_implementation
using ...Algorithms: CopyingBuffers

struct Input <: InputPart end

@kwdef mutable struct State <: StatePart
    sediment_buffer::Array{Float64, 4}
end

Common.initial_state(::Input, input::AbstractInput) = State(
    sediment_buffer=zeros(Float64, input.sediment_buffer_size, n_facies(input), input.box.grid_size...))

push_sediment!(s::State, args...) = CopyingBuffers.push_sediment!(s.sediment_buffer, args...)
pop_sediment!(s::State, args...) = CopyingBuffers.pop_sediment!(s.sediment_buffer, args...)
peek_sediment(s::State, args...) = CopyingBuffers.peek_sediment(s.sediment_buffer, args...)

default_implementation(::Type{AbstractSedimentBuffer}) = CopyingSedimentBuffer

end
```

To use the circular buffer, add `CircularSedimentBuffer` to the components of a model, for example `record_type(:Input, ALCAP, CircularSedimentBuffer)`. The initial state then follows the composition of the input.

``` {.julia file=src/Components/CircularSedimentBuffer.jl}
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
```

### Component Test

The component test runs for both implementations.

``` {.julia file=test/Components/SedimentBufferSpec.jl}
module SedimentBufferSpec
using Test
using Unitful
using CarboKitten
using CarboKitten.Components.Common: Amount, record_type
using CarboKitten.Components.Common: initial_state
using CarboKitten.SedimentBuffers: peek_sediment
import CarboKitten.Components.SedimentBuffer as SB
import CarboKitten.Components: CopyingSedimentBuffer, CircularSedimentBuffer

@testset "Components/SedimentBuffer" begin
    @testset "$(impl)" for impl in (CopyingSedimentBuffer, CircularSedimentBuffer)
        input = record_type(:Input, SB, impl)(
            box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
            time = TimeProperties(Δt=1.0u"yr", steps=10),
            facies = [record_type(:Facies, SB, impl)()],
            sediment_buffer_size = 10,
            depositional_resolution = 1.0u"m")
        state = initial_state(input)

        @test size(state.sediment_thickness) == (10, 1)
        @test state.sediment_buffer isa fieldtype(impl.State, :sediment_buffer)

        push! = SB.push_sediment(input)
        pop! = SB.pop_sediment(input)

        push!(state, reshape((1:10) .* 0.5u"m" |> collect, (1, 10, 1)))
        @test peek_sediment(state, 10.0)[1, :, 1] ≈ (1:10) .* 0.5
        @test peek_sediment(state, CartesianIndex(3, 1), 1.0) ≈ [1.0]
        # no initial topography, no subsidence
        @test state.sediment_thickness ≈ state.bathymetry

        buffer = zeros(Amount, 1, 10, 1)
        pop!(state, state.sediment_thickness ./ 2, buffer)
        @test reshape(state.sediment_thickness, (1, 10, 1)) ≈ buffer
        @test peek_sediment(state, 10.0)[1, :, 1] ≈ (1:10) .* 0.25
        # no initial topography, no subsidence
        @test state.sediment_thickness ≈ state.bathymetry
    end

    @testset "copying layout" begin
        input = record_type(:Input, SB)(
            box = Box{Periodic{2}}(grid_size=(10, 1), phys_scale=1.0u"m"),
            time = TimeProperties(Δt=1.0u"yr", steps=10),
            facies = [record_type(:Facies, SB)()],
            sediment_buffer_size = 10,
            depositional_resolution = 1.0u"m")
        state = initial_state(input)

        @test size(state.sediment_buffer) == (10, 1, 10, 1)

        SB.push_sediment(input)(state, reshape((1:10) .* 0.5u"m" |> collect, (1, 10, 1)))
        @test state.sediment_buffer[1, 1, :, 1] == repeat([0.5, 0.0], 5)

        SB.pop_sediment(input)(state, state.sediment_thickness ./ 2, zeros(Amount, 1, 10, 1))
        @test state.sediment_buffer[1, 1, :, 1] .% 1.0 ≈ repeat([0.25, 0.5, 0.75, 0.0], 3)[1:10]
    end
end
end
```
