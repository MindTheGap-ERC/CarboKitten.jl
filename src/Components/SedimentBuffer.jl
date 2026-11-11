# ~/~ begin <<docs/src/components/sediment_buffer.md#src/Components/SedimentBuffer.jl>>[init]
@component module SedimentBuffer
@mixin Boxes, FaciesBase, WaterDepth

using StaticArrays
using Unitful

using ..Common
using ..FaciesBase: n_facies
using CarboKitten.SedimentStack: pop_sediment!, push_sediment!, peek_sediment

export pop_sediment, push_sediment

@kwdef struct Input <: InputPart
    sediment_buffer_size::Int = 50
    depositional_resolution::Amount = 0.5u"m"
end

@kwdef mutable struct State <: StatePart
    sediment_buffer::Array{Float64,4}
    sediment_thickness::Array{Height, 2}
end

function push_sediment(input::AbstractInput)
    res = input.depositional_resolution
    n_f = n_facies(input)
    n_g = input.box.grid_size

    function (state::AbstractState, sediment::Array{Amount, 3})
        for i in CartesianIndices(n_g)
            total = sum(@view sediment[:, i[1], i[2]])
            state.sediment_thickness[i] += total
            state.bathymetry[i] += total
            v = SVector{n_f, Float64}(sediment[:, i[1], i[2]] ./ res .|> NoUnits)
            push_sediment!(view(state.sediment_buffer, :, :, i[1], i[2]), v)
        end
    end
end

function pop_sediment(input::AbstractInput)
    res = input.depositional_resolution
    n_g = input.box.grid_size

    function (state::AbstractState, amount::Array{Amount, 2}, out::Array{Amount, 3})
        for i in CartesianIndices(n_g)
            state.sediment_thickness[i] -= amount[i]
            state.bathymetry[i] -= amount[i]
            v = pop_sediment!(@view(state.sediment_buffer[:, :, i[1], i[2]]), amount[i] ./ res .|> NoUnits)
            view(out, :, i[1], i[2]) .= v .* res
        end
    end
end

end
# ~/~ end
