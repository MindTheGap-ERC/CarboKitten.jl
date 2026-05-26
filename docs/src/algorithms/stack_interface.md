# Sediment Stack Interface

We have two algorithms that provide a sediment buffer implementation.

``` {.julia file=src/Interfaces/SedimentBuffers.jl}
module SedimentBuffers

export push_sediment!, pop_sediment!, peek_sediment
using ..Chunks

"""
    push_sediment!(buffer, sediment)
    push_sediment!(buffer, sediment, chunk)

Push deposition onto the sediment buffer.
"""
function push_sediment! end

push_sediment!(buffer, sediment) = push_sediment!(buffer, sediment, Serial(Chunk(:, :)))

"""
    pop_sediment!(buffer, amount, out)
    pop_sediment!(buffer, amount, out, chunk)

Pop amount off the buffer, writing to sediment array.
"""
function pop_sediment! end

pop_sediment!(buffer, amount, out) = pop_sediment(buffer, amount, out, Serial(Chunk(:, :)))

"""
    peek_sediment(buffer, amount, out)
    peek_sediment(buffer, amount, out, chunk)

Peek at the top of the buffer.
"""
function peek_sediment end

peek_sediment(buffer, amount, out) = peek_sediment(buffer, amount, out, Serial(Chunk(:, :)))

end
```

## Tests

``` {.julia file=test/SedimentStackSpec.jl}
@testset "SedimentStack" begin
  using CarboKitten.SedimentStack: push_sediment!, pop_sediment!
  stack = zeros(Float64, 10, 3)
  @test pop_sediment!(stack, 0.0) == [0.0, 0.0, 0.0]
  push_sediment!(stack, [5.0, 0, 0])
  @test pop_sediment!(stack, 1.5) == [1.5, 0.0, 0.0]
  push_sediment!(stack, [0.0, 2.0, 0.0])   # (0 0.5) (0 1) (0.5 0.5) (1 0) ...
  @test pop_sediment!(stack, 2.0) == [0.25, 1.75, 0.0]
  @test pop_sediment!(stack, 1.5) == [1.25, 0.25, 0.0]
  @test pop_sediment!(stack, 0.0) == [0.0, 0.0, 0.0]
end

@testset "SedimentArray" begin
  using CarboKitten.SedimentStack: push_sediment!, peek_sediment
  sediment = zeros(Float64, 10, 3, 5, 5)
  for x in 1:10
    production = rand(3, 5, 5)
    push_sediment!(sediment, production)
  end
  a = peek_sediment(sediment, 1.0)
  @test all(sum(a; dims=1) .≈ 1.0)
end
```
