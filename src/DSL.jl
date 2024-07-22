# ~/~ begin <<docs/src/dsl.md#src/DSL.jl>>[init]
module DSL

include("DSL/Forward.jl")

using .Forward: @dynamic, @forward
using MacroTools: @capture, postwalk, prewalk
export @spec, @requires, @compose, @dynamic, @forward

# ~/~ begin <<docs/src/dsl.md#dsl-struct-type>>[init]
struct Struct
    mut::Bool
    kwarg::Bool
    parent::Union{Symbol, Nothing}
    fields::Vector{Union{Expr,Symbol}}
end

function define_struct(name::Symbol, s::Struct)
    if s.parent !== nothing
        name = :($name <: $(s.parent))
    end
    if s.mut
        :(mutable struct $name
            $(s.fields...)
        end)
    elseif s.kwarg
        :(@kwdef struct $name
            $(s.fields...)
          end)
    else
        :(struct $name
            $(s.fields...)
        end)
    end
end
# ~/~ end

# ~/~ begin <<docs/src/dsl.md#dsl>>[init]
"""
    @spec name body

Create a spec. When a spec is composed, the items in the spec will be spliced into a newly generated module. The `@spec` macro itself doesn't perform any operations other than storing the spec in a `const` expression. The real magic happens inside the `@compose` macro.
"""
macro spec(name, body)
    quoted_body = QuoteNode(body)

    clean_body = postwalk(e -> @capture(e, @requires parents__) ? :() : e, body)
    esc(Expr(:toplevel, :(module $name
        $(clean_body.args...)
        const AST = $quoted_body
    end)))
end

macro requires(deps...)
    esc(:(const PARENTS = [$(deps)...]))
end
# ~/~ end
# ~/~ begin <<docs/src/dsl.md#dsl>>[1]
macro compose(modname, cs, body)
    components = Set{Symbol}()

    structs = IdDict()
    using_statements = []
    const_statements = []
    specs_used = Set()

    # ~/~ begin <<docs/src/dsl.md#dsl-compose>>[init]
    function extend_struct!(name::Symbol, fields::Vector)
        append!(structs[name].fields, fields)
    end

    function create_struct!(name::Symbol, is_mutable::Bool, is_kwarg::Bool, abst::Union{Symbol, Nothing}, fields::Vector)
        structs[name] = Struct(is_mutable, is_kwarg, abst, fields)
    end

    function pass(e)
        if @capture(e, @requires parents__)
            parents .|> scan
            return
        end

        if @capture(e, (struct name_ fields__ end) |
                       (@kwdef struct kw_name_ fields__ end) |
                       (mutable struct mut_name_ fields__ end))
            is_mutable = mut_name !== nothing
            is_kwarg = kw_name !== nothing
            sname = is_mutable ? mut_name : (is_kwarg ? kw_name : name)

            @capture(sname, (name_ <: abst_) | name_)

            if name in keys(structs)
                extend_struct!(name, fields)
            else
                create_struct!(name, is_mutable, is_kwarg, abst, fields)
            end
            return
        end

        if @capture(e, const n_ = x_)
            push!(const_statements, e)
            return
        end

        if @capture(e, using x__ | using mod__: x__)
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

        e = Core.eval(__module__, :($(c).AST))
        prewalk(pass, e)
    end
    # ~/~ end

    @assert cs.head == :vect
    cs.args .|> scan

    Expr(:toplevel, esc(:(module $modname
        $(using_statements...)
        $(const_statements...)
        $(Iterators.map(splat(define_struct), pairs(structs))...)
        $(body.args...)
    end)))
end
# ~/~ end

end
# ~/~ end
