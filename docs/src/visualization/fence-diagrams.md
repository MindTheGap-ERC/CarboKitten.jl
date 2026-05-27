# Fence Diagrams

The fence-diagram style visualization allows users to generate multiple cross-sections through the model domain by defining the number of fences to display, along with their x and y positions. This provides an intuitive way to inspect the 3D stratigraphic architecture and examine both lateral and vertical variations in facies, sediment accumulation, or other stratigraphic properties across the platform.

``` {.julia file="ext/FenceDiagrams.jl"}
module FenceDiagrams

import CarboKitten.Visualization: fence_diagrams, fence_diagrams!

end
```
