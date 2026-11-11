# ~/~ begin <<docs/src/algorithms/circular_stack.md#test/Algorithms/CircularBuffersSpec.jl>>[init]
module CircularBuffersSpec
using Test
using CarboKitten.Algorithms.CircularBuffers: CircularBuffer, push_sediment!, pop_sediment!, peek_sediment, unroll
using CarboKitten.Algorithms: CopyingBuffers
using Random: Xoshiro

@testset "Algorithms/CircularBuffers" begin
  @testset "column" begin
    col = zeros(Float64, 10, 3)
    head = 1
    parcel = zeros(Float64, 3)
    head = pop_sediment!(col, head, 0.0, parcel)
    @test parcel == [0.0, 0.0, 0.0]
    head = push_sediment!(col, head, [5.0, 0, 0])
    head = pop_sediment!(col, head, 1.5, parcel)
    @test parcel == [1.5, 0.0, 0.0]
    head = push_sediment!(col, head, [0.0, 2.0, 0.0])
    head = pop_sediment!(col, head, 2.0, parcel)
    @test parcel == [0.25, 1.75, 0.0]
    head = pop_sediment!(col, head, 1.5, parcel)
    @test parcel == [1.25, 0.25, 0.0]
    head = pop_sediment!(col, head, 0.0, parcel)
    @test parcel == [0.0, 0.0, 0.0]
  end

  @testset "edge cases" begin
    # pushing more than the buffer can hold fills every row with the parcel's fractions
    col = zeros(Float64, 4, 2)
    head = @test_logs (:warn,) push_sediment!(col, 1, [6.0, 2.0])
    @test all(col .== [0.75 0.25])

    # popping more than the buffer holds returns everything
    col = [0.0 1.0; 0.0 1.0; 1.0 0.0; 1.0 0.0]
    parcel = zeros(Float64, 2)
    @test_logs (:error,) pop_sediment!(col, 4, 10.0, parcel)
    @test parcel == [2.0, 2.0]
    @test all(col .== 0.0)

    # the head wraps around
    col = zeros(Float64, 4, 2)
    head = push_sediment!(col, 4, [1.5, 0.0])
    @test head == 1
    @test col[:, 1] == [0.5, 0.0, 0.0, 1.0]
    @test unroll(col, head) == [0.5 0.0; 1.0 0.0; 0.0 0.0; 0.0 0.0]
  end

  @testset "equivalent to copying stack" begin
    rng = Xoshiro(0)
    for _ in 1:20
      copying = zeros(Float64, 8, 3)
      circular = zeros(Float64, 8, 3)
      head = rand(rng, 1:8)
      for _ in 1:50
        if rand(rng, Bool)
          parcel = rand(rng, 3) .* rand(rng, (0.1, 1.0, 2.0))
          CopyingBuffers.push_sediment!(copying, parcel)
          head = push_sediment!(circular, head, parcel)
        else
          Δ = rand(rng) * rand(rng, (0.1, 1.0, 3.0))
          @test peek_sediment(circular, head, Δ) ≈ CopyingBuffers.peek_sediment(copying, Δ)
          expected = CopyingBuffers.pop_sediment!(copying, Δ)
          parcel = zeros(Float64, 3)
          head = pop_sediment!(circular, head, Δ, parcel)
          @test parcel ≈ expected
        end
        @test unroll(circular, head) ≈ copying
      end
    end
  end

  @testset "grid" begin
    buffer = CircularBuffer{Float64}(10, 3, (5, 5))
    for x in 1:10
      production = rand(3, 5, 5)
      push_sediment!(buffer, production)
    end
    a = peek_sediment(buffer, 1.0)
    @test all(sum(a; dims=1) .≈ 1.0)

    out = zeros(Float64, 3, 5, 5)
    pop_sediment!(buffer, fill(1.0, 5, 5), out)
    @test out ≈ a
  end
end
end
# ~/~ end
