# Box

This module makes sure we have access to box properties.

``` {.julia file=src/Components/Boxes.jl}
@compose module Boxes
    using ..Common

    @kwdef struct Input <: AbstractInput
        box::Box
    end
end
```
