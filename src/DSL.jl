# ~/~ begin <<docs/src/dsl.md#src/DSL.jl>>[init]
module DSL

using MacroTools: @capture, postwalk

# ~/~ begin <<docs/src/dsl.md#dsl>>[init]
"""
    @spec name body

Create a spec. When a spec is composed, the items in the spec will be spliced into a newly generated module. The `@spec` macro itself doesn't perform any operations other than storing the spec in a `const` expression.
"""
macro spec(name, body)
	:(const $(esc(name)) = $(QuoteNode(body)))
end
# ~/~ end
# ~/~ begin <<docs/src/dsl.md#dsl>>[1]
# ~/~ begin <<docs/src/dsl.md#dsl-struct-type>>[init]
struct Struct
	mut::Bool
	fields::Vector{Union{Expr,Symbol}}
end

function define_struct(name::Symbol, s::Struct)
	if s.mut
		:(mutable struct $name
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

  function create!(name::Symbol, is_mutable::Bool, fields::Vector)
  	structs[name] = Struct(is_mutable, fields)
  end

  function pass(e)
  	if @capture(e, @requires parents__)
  		parents .|> scan
  		return e
  	end

  	if @capture(e, (struct name_ fields__ end) |
  				   (mutable struct mut_name_ fields__ end))
  		is_mutable = mut_name !== nothing
  		name = is_mutable ? mut_name : name

  		if name in keys(structs)
  			extend!(name, fields)
  		else
  			create!(name, is_mutable, fields)
  		end
  		return e
  	end

  	if @capture(e, const n_ = x_)
  		const_statements[n] = x
  		return e
  	end

  	if @capture(e, using x__)
  		push!(using_statements, e)
  		return e
  	end
  end

  function scan(c::Symbol)
  	if c in specs_used
  		return
  	end
  	push!(specs_used, c)

  	postwalk(pass, esc(c))
  end
  # ~/~ end

  @assert cs.head == :vect
  cs.args .|> scan

  :(module $(esc(modname))
	  $(using_statements...)
	  $(Iterators.map(splat(define_const), pairs(const_statements))...)
	  $(Iterators.map(splat(define_struct), pairs(structs))...)
	end)
end
# ~/~ end

end
# ~/~ end
