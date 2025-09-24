# Tags

```component-dag
CarboKitten.Components.Tag
```

This module allows simulations to be tagged with a string.

``` {.julia file=src/Components/Tag.jl}
@compose module Tag
using ..Common
using HDF5

@kwdef struct Input <: AbstractInput
    tag::String = "untagged run"
end

function write_header(input::AbstractInput, output::AbstractOutput)
    set_attribute(output, "tag", input.tag)
end
end
```
