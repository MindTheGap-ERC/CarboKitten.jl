# ~/~ begin <<docs/src/components/sediment_buffer.md#test/SedimentStackSpec.jl>>[init]
@testset "SedimentStack" begin
  using CarboKitten.SedimentStack: push_sediment!, pop_sediment!, peek_sediment
  stack = zeros(Float64, 10, 3)
  @test pop_sediment!(stack, 0.0) == [0.0, 0.0, 0.0]
  push_sediment!(stack, [5.0, 0, 0])
  @test pop_sediment!(stack, 1.5) == [1.5, 0.0, 0.0]
  push_sediment!(stack, [0.0, 2.0, 0.0])   # (0 0.5) (0 1) (0.5 0.5) (1 0) ...
  @test pop_sediment!(stack, 2.0) == [0.25, 1.75, 0.0]
  @test pop_sediment!(stack, 1.5) == [1.25, 0.25, 0.0]
  @test pop_sediment!(stack, 0.0) == [0.0, 0.0, 0.0]

  stack2 = zeros(Float64, 10, 3)
  push_sediment!(stack2, [1.0, 0.0, 0.0])
  @test peek_sediment(stack2, 1.0) == [1.0, 0.0, 0.0]
end

@testset "SedimentStack edge cases" begin
  using CarboKitten.SedimentStack: push_sediment!, pop_sediment!

  # pushing more than the buffer can hold fills every row with the parcel's fractions
  stack = zeros(Float64, 4, 2)
  @test_logs (:warn,) push_sediment!(stack, [6.0, 2.0])
  @test all(stack .== [0.75 0.25])

  # popping several rows from a full buffer keeps the deepest rows
  stack = [1.0 0.0; 1.0 0.0; 0.0 1.0; 0.0 1.0]
  @test pop_sediment!(stack, 2.0) == [2.0, 0.0]
  @test stack == [0.0 1.0; 0.0 1.0; 0.0 0.0; 0.0 0.0]

  # popping more than the buffer holds returns everything, counting each row once
  stack = [1.0 0.0; 1.0 0.0; 0.0 1.0; 0.0 1.0]
  @test (@test_logs (:error,) pop_sediment!(stack, 10.0)) == [2.0, 2.0]
  @test all(stack .== 0.0)
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
# ~/~ end
