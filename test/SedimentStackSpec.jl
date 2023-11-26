# ~/~ begin <<docs/src/submarine-transport.md#test/SedimentStackSpec.jl>>[init]
@testset "SedimentStack" begin
  stack = zeros(Float64, 10, 3)
  push_sediment!(stack, [5.0, 0, 0])
  @test pop_sediment!(stack, 1.5) == [1.5, 0.0, 0.0]
end
# ~/~ end