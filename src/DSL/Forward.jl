# ~/~ begin <<docs/src/dsl.md#src/DSL/Forward.jl>>[init]
module Forward

using MacroTools: @capture, postwalk, prewalk
import Base: getproperty, setproperty!

struct Field{Sym} end

function define_dynamic(class::Symbol)
    getter = Symbol(lowercase(string(class)), "_get")
    setter = Symbol(lowercase(string(class)), "_set!")
    :(begin
        $(esc(getter))(self::$(esc(class)), ::$(esc(:Type)){Field{symb}}) where {symb} =
            $(esc(:(Core.getfield)))(self, symb)
        $(esc(setter))(self::$(esc(class)), ::$(esc(:Type)){Field{symb}}, value) where {symb} =
            $(esc(:(Core.setfield!)))(self, symb, value)
        $(esc(:(Base.getproperty)))(self::$(esc(class)), symb::Symbol) =
            $(esc(getter))(self, Field{symb})
        $(esc(:(Base.setproperty!)))(self::$(esc(class)), symb::Symbol, value) =
            $(esc(setter))(self, Field{symb}, value)
    end)
end

function get_path_expr(path...)
    foldl((s, f) -> :($(esc(:(Core.getfield)))($s, $(Expr(:quote, f)))), path)
end

function define_forward(class::Symbol, name::Symbol, path::Tuple)
    getter = Symbol(lowercase(string(class)), "_get")
    setter = Symbol(lowercase(string(class)), "_set!")
    :(begin
        $(esc(getter))(self::$(esc(class)), ::$(esc(:Type)){Field{$(Expr(:quote, name))}}) =
            $(get_path_expr(:self, path...))
        $(esc(setter))(self::$(esc(class)), ::$(esc(:Type)){Field{$(Expr(:quote, name))}}, value) =
            $(esc(:(Core.setfield!)))($(get_path_expr(:self, Base.front(path)...)), $(Expr(:quote, last(path))), value)
    end)
end

macro dynamic(class)
    define_dynamic(class)
end

function get_path(expr)
    if @capture(expr, a_.b_)
        return (get_path(a)..., b)
    else
        return (expr,)
    end
end

macro forward(expr)
    @assert @capture(expr, fwd_ ~ path_)
    @assert @capture(fwd, class_.field_)
    define_forward(class, field, Base.tail(get_path(path)))
end

end
# ~/~ end
