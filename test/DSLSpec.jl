# ~/~ begin <<docs/src/dsl.md#test/DSLSpec.jl>>[init]
using CarboKitten.DSL
using MacroTools: prewalk, rmlines

clean(expr) = prewalk(rmlines, expr)

# ~/~ begin <<docs/src/dsl.md#dsl-spec-defs>>[init]
struct Parent
  x::Int
end

struct Child
  p::Parent
  y::Int
end

@dynamic Child
@forward Child.x ~ Child.p.x
# ~/~ end
# ~/~ begin <<docs/src/dsl.md#dsl-spec-defs>>[1]
@spec MySpec begin
    const msg = "hello"
end
# ~/~ end
# ~/~ begin <<docs/src/dsl.md#dsl-spec-defs>>[2]
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
# ~/~ begin <<docs/src/dsl.md#dsl-spec-defs>>[3]
@spec C begin
  @requires A
  struct S
    c::Int
  end

  @kwdef struct T
    f::Int
  end
end

@compose AC [C]
# ~/~ end

@testset "CarboKitten.DSL" begin
  # ~/~ begin <<docs/src/dsl.md#dsl-spec>>[init]
  let c = Child(Parent(42), 23)
  @assert c.x == 42
  @assert c.y == 23
  end
  # ~/~ end
  # ~/~ begin <<docs/src/dsl.md#dsl-spec>>[1]
  @test clean(MySpec.AST) == clean(:(begin const msg = "hello" end))
  @test MySpec.msg == "hello"
  # ~/~ end
  # ~/~ begin <<docs/src/dsl.md#dsl-spec>>[2]
  @test fieldnames(AB.S) == (:a, :b)
  # ~/~ end
  # ~/~ begin <<docs/src/dsl.md#dsl-spec>>[3]
  @test fieldnames(AC.S) == (:a, :c)
  @test fieldnames(AC.T) == (:f,)
  @test AC.T(f = 4).f == 4
  # ~/~ end
end
# ~/~ end
