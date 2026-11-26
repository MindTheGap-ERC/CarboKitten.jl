# Custom Algorithms

``` {.julia file=test/AlgorithmsSpec.jl}
module AlgorithmsSpec
using Test
using CarboKitten.Algorithms: enumerate_seq

@testset "CarboKitten.Algorithms" begin
    <<algorithms-spec>>
end
end
```
