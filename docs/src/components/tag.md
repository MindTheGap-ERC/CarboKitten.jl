# Tags

```component-dag
CarboKitten.Components.Tag
```

This module allows simulations to be tagged with a string.

``` {.julia file=src/Components/Tag.jl}
@component module Tag
using ..Common
import ..Common: write_header
using HDF5

@kwdef struct Input <: InputPart
    tag::String = "untagged run"
end

function write_header(::Input, input::AbstractInput, output::AbstractOutput)
    set_attribute(output, "tag", input.tag)
end
end
```
