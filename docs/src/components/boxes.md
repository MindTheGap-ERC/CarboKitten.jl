``` {.julia file=src/Components/Boxes.jl}
@compose module Boxes
    using ..Common

    @kwdef struct Input <: AbstractInput
        box::Box
    end
end
```
