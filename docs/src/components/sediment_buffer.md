# Sediment Buffer

``` {.julia file=src/Components/SedimentBuffer.jl}
@compose module SedimentBuffer
@mixin Boxes

using ..Common
using CarboKitten.SedimentStack: pop_sediment!, push_sediment!, peek_sediment

export pop_sediment!, push_sediment!, peek_sediment

@kwdef struct Input <: AbstractInput
    sediment_buffer_size::Int = 50
    depositional_resolution::Amount = 0.5u"m"
end

@kwdef mutable struct State <: AbstractState
    sediment_buffer::Array{Float64,4}
end

end
```

