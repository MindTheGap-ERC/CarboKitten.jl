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

function write_header(fid, input::AbstractInput)
    attr = attributes(fid["input"])
    attr["tag"] = input.tag
end
end
```
