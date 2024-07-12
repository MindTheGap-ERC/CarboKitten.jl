# The CarboKitten DSL

At some point in the evolution of CarboKitten, there will be many different implementations for production, disintegration, transport etc. Each with slightly different needs on the `Input` and `State` structures. Wouldn't it be cool if we could compose a model directly out of these components?

```@raw html
<details><summary>Boiler plate</summary>
```

``` {.julia file=test/DSLSpec.jl}
using CarboKitten.DSL: @spec, @compose
using MacroTools: prewalk, rmlines

clean(expr) = prewalk(rmlines, expr)

<<dsl-spec-defs>>

@testset "CarboKitten.DSL" begin
  <<dsl-spec>>
end
```

``` {.julia file=src/DSL.jl}
module DSL

using MacroTools: @capture, postwalk

<<dsl>>

end
```

```@raw html
</details>
```

## Example

As a first example, let us recreate the [Bosscher1992](@cite) model. For this, we need to know the water depth, and specify a uniform production (i.e. without CA). First we define some units.

``` {.julia file=examples/dsl/bs92.jl}
using CarboKitten.DSL: @spec, @compose

module Units
  using Unitful
  using CarboKitten.Config: Box

  export Amount, Time, Height, Location, Rate, Intensity, Box

  const Amount = typeof(1.0u"m")
  const Time = typeof(1.0u"Myr")
  const Height = typeof(1.0u"m")
  const Location = typeof(1.0u"km")
  const Rate = typeof(1.0u"m/Myr")
  const Intensity = typeof(1.0u"W/m^2")
end

<<dsl-example-time>>
<<dsl-example-waterdepth>>
<<dsl-example-production>>

@compose BS92 [UniformProduction]

module BS92

end
```

### Time

``` {.julia #dsl-example-time}
@spec TimeIntegration begin
  using ..Types
  using CarboKitten.Config: TimeProperties

  struct Input
    time::TimeProperties
  end

  struct State
    time::Time
  end
end
```

### Water depth
``` {.julia #dsl-example-waterdepth}
@spec WaterDepth begin
  @requires TimeIntegration
  using ..Types

  struct Input
    box::Box
    sea_level          # function (t::Time) -> Length
    bedrock_elevation  # function (x::Location, y::Location) -> Length
    subsidence_rate::Rate
  end

  struct State
    time::Time
    sediment_height::Matrix{Height}
  end
end

function water_depth(input)
  x, y = axes(input.box)
  eta0 = input.bedrock_elevation.(x, y')

  return function(state::State)
    return input.sea_level(state.time) .- eta0 .+
      (input.subsidence_rate * state.time) .- state.sediment_height
  end
end
```

### Uniform Production

``` {.julia #dsl-example-production}
@spec UniformProduction begin
  @requires WaterDepth
  using ..Types

  struct Facies
    maximum_growth_rate::Rate
    extinction_coefficient::typeof(1.0u"m^-1")
    saturation_intensity::Intensity
  end

  struct Input
    insolation::Intensity
    facies::Vector{Facies}
  end
end

function production_rate(insolation, facies, water_depth)
    gₘ = facies.maximum_growth_rate
    I = insolation / facies.saturation_intensity
    x = water_depth * facies.extinction_coefficient
    return water_depth > 0.0u"m" ? gₘ * tanh(I * exp(-x)) : 0.0u"m/Myr"
end

function uniform_production(input)
  w = water_depth(input)
  na = [CartesianIndex()]

  return function(state)
    return production_rate.(
      input.insolation,
      input.facies[:,na,na],
      w(state)[na,:,:])
  end
end
```
## `@spec`

The `@spec` macro stores a spec syntax.

``` {.julia #dsl-spec-defs}
@spec MySpec begin
    "hello"
end
```

``` {.julia #dsl-spec}
@test clean(MySpec) == clean(:(begin "hello" end))
```

The `@spec` macro is used to specify the structs of a model component.

``` {.julia #dsl}
"""
    @spec name body

Create a spec. When a spec is composed, the items in the spec will be spliced into a newly generated module. The `@spec` macro itself doesn't perform any operations other than storing the spec in a `const` expression. The real magic happens inside the `@compose` macro.
"""
macro spec(name, body)
    :(const $(esc(name)) = $(QuoteNode(body)))
end
```

## `@compose`

The idea of `@compose` is that it splices `struct` definitions, such that resulting structs contain all members from required specs.

We define some variables to collect structs, consts and `using` declarations. At the end we use these collections to build a new module.

``` {.julia #dsl-spec-defs}
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
```

``` {.julia #dsl-spec}
@test fieldnames(AB.S) == (:a, :b)
```

A spec can depend on another using the `@require` syntax.

``` {.julia #dsl-spec-defs}
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
```

``` {.julia #dsl-spec}
@test fieldnames(AC.S) == (:a, :c)
@test fieldnames(AC.T) == (:f,)
@test AC.T(f = 4).f == 4
```

```@raw html
<details><summary>`@compose` implementation</summary>
```

``` {.julia #dsl}
<<dsl-struct-type>>

function define_const(name::Symbol, v)
    :(const $(esc(name)) = $v)
end

macro compose(modname, cs)
    components = Set{Symbol}()

    structs = IdDict()
    using_statements = []
    const_statements = IdDict()
    specs_used = Set()

    <<dsl-compose>>

    @assert cs.head == :vect
    cs.args .|> scan

    Core.eval(__module__, :(module $modname
        $(using_statements...)
        $(Iterators.map(splat(define_const), pairs(const_statements))...)
        $(Iterators.map(splat(define_struct), pairs(structs))...)
    end))
end
```

``` {.julia #dsl-compose}
function extend!(name::Symbol, fields::Vector)
    append!(structs[name].fields, fields)
end

function create!(name::Symbol, is_mutable::Bool, is_kwarg::Bool, fields::Vector)
    structs[name] = Struct(is_mutable, is_kwarg, fields)
end

function pass(e)
    if @capture(e, @requires parents__)
        parents .|> scan
        return
    end

    if @capture(e, (struct name_ fields__ end) |
                   (@kwdef struct kw_name_ fields__ end)
                   (mutable struct mut_name_ fields__ end))
        is_mutable = mut_name !== nothing
        is_kwarg = kw_name !== nothing
        name = is_mutable ? mut_name : (is_kwarg ? kw_name : name)

        if name in keys(structs)
            extend!(name, fields)
        else
            create!(name, is_mutable, is_kwarg, fields)
        end
        return
    end

    if @capture(e, const n_ = x_)
        const_statements[n] = x
        return
    end

    if @capture(e, using x__)
        push!(using_statements, e)
        return
    end

    return e
end

function scan(c::Symbol)
    if c in specs_used
        return
    end
    push!(specs_used, c)

    e = Core.eval(__module__, :($c))
    postwalk(pass, e)
end
```

``` {.julia #dsl-struct-type}
struct Struct
    mut::Bool
    kwarg::Bool
    fields::Vector{Union{Expr,Symbol}}
end

function define_struct(name::Symbol, s::Struct)
    if s.mut
        :(mutable struct $name
            $(s.fields...)
        end)
    elseif s.kwarg
        :(@kwarg struct $name
            $(s.fields...)
          end)
    else
        :(struct $name
            $(s.fields...)
        end)
    end
end
```

```@raw html
</details>
```

## Some types

``` {.julia #types}

abstract type Input end
abstract type State end
abstract type Facies end
```

