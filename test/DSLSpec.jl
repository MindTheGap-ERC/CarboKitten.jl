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

@compose AB [A, B]
# ~/~ end
# ~/~ begin <<docs/src/dsl.md#dsl-spec-defs>>[2]
@spec C begin
  @requires A
  struct S
    c::Int
  end

  @kwarg struct T
    f::Int
  end
end

@compose AC [C]
# ~/~ end

@testset "CarboKitten.DSL" begin
  # ~/~ begin <<docs/src/dsl.md#dsl-spec>>[init]
  @test clean(MySpec) == clean(:(begin "hello" end))
  # ~/~ end
  # ~/~ begin <<docs/src/dsl.md#dsl-spec>>[1]
  @test fieldnames(AB.S) == (:a, :b)
  # ~/~ end
  # ~/~ begin <<docs/src/dsl.md#dsl-spec>>[2]
  @test fieldnames(AC.S) == (:a, :c)
  @test fieldnames(AC.T) == (:f,)
  @test AC.T(f = 4).f == 4
  # ~/~ end
end
# ~/~ end
