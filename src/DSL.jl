# ~/~ begin <<docs/src/dsl.md#src/DSL.jl>>[init]
module DSL

using MacroTools: @capture, postwalk

# ~/~ begin <<docs/src/dsl.md#dsl>>[init]
"""
    @spec name body

Create a spec. When a spec is composed, the items in the spec will be spliced into a newly generated module. The `@spec` macro itself doesn't perform any operations other than storing the spec in a `const` expression. The real magic happens inside the `@compose` macro.
"""
macro spec(name, body)
    :(const $(esc(name)) = $(QuoteNode(body)))
end
# ~/~ end
# ~/~ begin <<docs/src/dsl.md#dsl>>[1]
# ~/~ begin <<docs/src/dsl.md#dsl-struct-type>>[init]
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
# ~/~ end

function define_const(name::Symbol, v)
    :(const $(esc(name)) = $v)
end

macro compose(modname, cs)
    components = Set{Symbol}()

    structs = IdDict()
    using_statements = []
    const_statements = IdDict()
    specs_used = Set()

    # ~/~ begin <<docs/src/dsl.md#dsl-compose>>[init]
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
    # ~/~ end

    @assert cs.head == :vect
    cs.args .|> scan

    Core.eval(__module__, :(module $modname
        $(using_statements...)
        $(Iterators.map(splat(define_const), pairs(const_statements))...)
        $(Iterators.map(splat(define_struct), pairs(structs))...)
    end))
end
# ~/~ end

end
# ~/~ end
