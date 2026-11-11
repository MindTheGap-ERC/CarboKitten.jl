# ~/~ begin <<docs/src/components/sediment_buffer.md#src/Components/SedimentBuffer.jl>>[init]
@component module SedimentBuffer
@mixin Boxes, FaciesBase, WaterDepth
@requires AbstractSedimentBuffer

using Unitful

using ..Common
using ..FaciesBase: n_facies
using ...SedimentBuffers: AbstractSedimentBuffer, initial_sediment_buffer, push_sediment!, pop_sediment!

export pop_sediment, push_sediment, initial_sediment_buffer

@kwdef struct Input <: InputPart
    sediment_buffer_size::Int = 50
    depositional_resolution::Amount = 0.5u"m"
end

@kwdef mutable struct State <: StatePart
    sediment_thickness::Array{Height, 2}
end

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
# ~/~ end
