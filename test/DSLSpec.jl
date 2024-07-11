# ~/~ begin <<docs/src/dsl.md#test/DSLSpec.jl>>[init]
using CarboKitten.DSL: @spec, @compose
using MacroTools: prewalk, rmlines

clean(expr) = prewalk(rmlines, expr)

# ~/~ begin <<docs/src/dsl.md#dsl-spec-defs>>[init]
@spec MySpec begin
    "hello"
end
# ~/~ end
# ~/~ begin <<docs/src/dsl.md#dsl-spec-defs>>[1]
@spec A begin
  struct S
    a::Int
  end
end

@spec B begin
  struct S
    b::Int
  end
end

@compose C [A, B]
# ~/~ end

@testset "CarboKitten.DSL" begin
  # ~/~ begin <<docs/src/dsl.md#dsl-spec>>[init]
  @test clean(MySpec) == clean(:(begin "hello" end))
  # ~/~ end
  # ~/~ begin <<docs/src/dsl.md#dsl-spec>>[1]
  @test getfields(C.S) == [:a, :b]
  # ~/~ end
end
# ~/~ end