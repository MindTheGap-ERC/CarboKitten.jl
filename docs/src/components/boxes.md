# Box

```component-dag
CarboKitten.Components.Boxes
```

This module makes sure we have access to box properties.

``` {.julia file=test/Components/BoxesSpec.jl}
module BoxesSpec
    using Test
    using CarboKitten.Components.Common
    using CarboKitten.Components.Boxes

    @testset "Components/Boxes" begin
        let box = Box{Periodic{2}}(grid_size=(10, 10), phys_scale=2.0u"m"),
            input = Boxes.Input(box=box)
            @test input.box.grid_size == (10, 10)
            @test input.box.phys_scale == 2.0u"m"
        end
    end
end
```

``` {.julia file=src/Components/Boxes.jl}
@compose module Boxes
using ..Common

@kwdef struct Input <: AbstractInput
    box::Box
end

function write_header(fid, input::AbstractInput)
    x, y = Common.axes(input.box)

    gid = fid["input"]
    gid["x"] = collect(x) |> in_units_of(u"m")
    gid["y"] = collect(y) |> in_units_of(u"m")
end
end
```
