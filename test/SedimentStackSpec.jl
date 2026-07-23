# ~/~ begin <<docs/src/components/sediment_buffer.md#test/SedimentStackSpec.jl>>[init]
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
@testset "Sediment layer reconstruction" begin
    using CarboKitten.SedimentStack: sediment_layer

    deposition = zeros(Float64, 3, 1, 1, 4)
    disintegration = zeros(Float64, 3, 1, 1, 4)

    deposition[1, 1, 1, 1] = 1.0
    deposition[2, 1, 1, 2] = 1.0
    deposition[3, 1, 1, 3] = 1.0

    top, present = sediment_layer(deposition, disintegration, 3)
    @test present[1, 1]
    @test top[:, 1, 1] ≈ [0.0, 0.0, 1.0]

    middle, present = sediment_layer(deposition, disintegration, 3; depth=1.0)
    @test present[1, 1]
    @test middle[:, 1, 1] ≈ [0.0, 1.0, 0.0]

    bottom, present = sediment_layer(deposition, disintegration, 3; depth=2.0)
    @test present[1, 1]
    @test bottom[:, 1, 1] ≈ [1.0, 0.0, 0.0]

    combined, present = sediment_layer(
        deposition,
        disintegration,
        3;
        thickness=2.0,
    )
    @test present[1, 1]
    @test combined[:, 1, 1] ≈ [0.0, 1.0, 1.0]

    earlier, present = sediment_layer(deposition, disintegration, 2)
    @test present[1, 1]
    @test earlier[:, 1, 1] ≈ [0.0, 1.0, 0.0]

    absent, present = sediment_layer(deposition, disintegration, 3; depth=3.0)
    @test !present[1, 1]
    @test iszero(sum(absent[:, 1, 1]))

    eroded_deposition = zeros(Float64, 2, 1, 1, 3)
    eroded_disintegration = zeros(Float64, 2, 1, 1, 3)
    eroded_deposition[1, 1, 1, 1] = 1.0
    eroded_deposition[2, 1, 1, 2] = 1.0
    eroded_disintegration[2, 1, 1, 3] = 1.0

    exposed, present = sediment_layer(
        eroded_deposition,
        eroded_disintegration,
        3,
    )
    @test present[1, 1]
    @test exposed[:, 1, 1] ≈ [1.0, 0.0]

    # A finite layer thickness combines a very small final event with the
    # underlying sediment instead of letting that event dominate the map.
    thin_deposition = zeros(Float64, 2, 1, 1, 2)
    thin_disintegration = zeros(Float64, 2, 1, 1, 2)
    thin_deposition[1, 1, 1, 1] = 1.0
    thin_deposition[2, 1, 1, 2] = 0.01

    smoothed, present = sediment_layer(
        thin_deposition,
        thin_disintegration,
        2;
        thickness=1.0,
    )
    @test present[1, 1]
    @test smoothed[:, 1, 1] ≈ [0.99, 0.01]

    # Coverage must be nested with depth: once a preserved column is too thin
    # at one depth, it must remain absent at every greater depth.
    shallow, present_shallow = sediment_layer(
        deposition,
        disintegration,
        3;
        depth=3.0,
    )
    deep, present_deep = sediment_layer(
        deposition,
        disintegration,
        3;
        depth=5.0,
    )
    @test !present_shallow[1, 1]
    @test !present_deep[1, 1]
    @test iszero(sum(shallow))
    @test iszero(sum(deep))
end
# ~/~ end
