### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ a5f0eccc-b53d-4d41-88cd-458ec07545ab
using Pkg; Pkg.activate("../workenv")

# ╔═╡ 11b1526c-a419-445b-8b74-f5b15fba5550
using Unitful

# ╔═╡ 18d7285e-3bc5-44c1-8ba6-8dc1becefa3d
using MacroTools: @capture, postwalk

# ╔═╡ f44c7bd5-a233-4605-bc1c-fe2a085cdbef
using CarboKitten.Config: Box, axes

# ╔═╡ d8b3f8a2-9110-4be0-8b86-a1d67610efe4
module Units
	using Unitful
	using CarboKitten.Config: Box

	export Amount, Time, Height, Location, Rate, Intensity, AbstractInput, AbstractState, AbstractFacies, Box

	const Amount = typeof(1.0u"m")
	const Time = typeof(1.0u"Myr")
	const Height = typeof(1.0u"m")
	const Location = typeof(1.0u"km")
	const Rate = typeof(1.0u"m/Myr")
	const Intensity = typeof(1.0u"W/m^2")
	
	abstract type AbstractInput end
	abstract type AbstractState end
	abstract type AbstractFacies end
end

# ╔═╡ 72e69f4a-9835-4f6b-b61d-0e5ff12319dd
macro compose(modname, cs...)
	structs = IdDict()
	used = Set{Symbol}()
	using_statements = []
	const_statements = IdDict()

	function extend!(name::Symbol, fields::Vector)
		append!(structs[name].fields, fields)
	end

	function create!(name::Symbol, is_mutable::Bool, fields::Vector)
		structs[name] = Struct(is_mutable, fields)
	end
	
	function pass(e)
		if @capture(e, @requires parents__)
			parents .|> use
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

	function use(c::Symbol)
		if c in used
			return
		end
		push!(used, c)
		
		postwalk(pass, eval(c))
	end

	function define_const(name::Symbol, v)
		:(const $name = $v)
	end
	
	function define(name::Symbol, s::Struct)
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
	
	cs .|> use

	@eval module $modname
		using Unitful
		
		$(using_statements...)
		$(Iterators.map(splat(define_const), pairs(const_statements)))
	   	$(Iterators.map(splat(define), pairs(structs))...)
	end
end

# ╔═╡ a644565f-fa1e-4040-817d-90ab36eacec9
macro spec(block)
	return QuoteNode(block)
end

# ╔═╡ 3777e60a-3e27-11ef-2e27-bf4856556e6d
const WaterDepth = @spec begin
  using CarboKitten.Config: Box
	const Time = typeof(1.0u"Myr")
	const Height = typeof(1.0u"m")
	const Rate = typeof(1.0u"m/Myr")
  
  struct Input
    box::Box
    sea_level          # function (t::Time) -> Length
    bedrock_elevation  # function (x::Location, y::Location) -> Length
    subsidence_rate::Rate
  end

  mutable struct State
    time::Time
    sediment_height::Matrix{Height}
  end
end;

# ╔═╡ bad996f3-6a72-497a-8437-b2e2259f284c
const UniformProduction = @spec begin
  @requires WaterDepth

	const Rate = typeof(1.0u"m/Myr")
	const Intensity = typeof(1.0u"W/m^2")
	
  struct Facies
    maximum_growth_rate::Rate
    extinction_coefficient::typeof(1.0u"m^-1")
    saturation_intensity::Intensity
  end

  struct Input
    insolation::Intensity
    facies::Vector{Facies}
  end
end;

# ╔═╡ 374ce450-dc6d-4b4a-8697-9926111686fc
struct Struct
	mut::Bool
	fields::Vector{Union{Expr,Symbol}}
end

# ╔═╡ aaa31f40-a910-4f93-9bd9-601a8814a449
UniformProduction.args[2].head

# ╔═╡ 8227fd27-33de-44e7-82d2-5065e5d368fb
module Test
end

# ╔═╡ a08ade9f-134f-43ed-8b32-17f7750c87b2
@compose UPMP UniformProduction

# ╔═╡ 5d793a94-7f96-4511-841a-de622dc50662
fieldnames(UPMP.Input)

# ╔═╡ Cell order:
# ╠═a5f0eccc-b53d-4d41-88cd-458ec07545ab
# ╠═11b1526c-a419-445b-8b74-f5b15fba5550
# ╠═18d7285e-3bc5-44c1-8ba6-8dc1becefa3d
# ╠═f44c7bd5-a233-4605-bc1c-fe2a085cdbef
# ╠═d8b3f8a2-9110-4be0-8b86-a1d67610efe4
# ╠═a644565f-fa1e-4040-817d-90ab36eacec9
# ╠═3777e60a-3e27-11ef-2e27-bf4856556e6d
# ╠═bad996f3-6a72-497a-8437-b2e2259f284c
# ╠═374ce450-dc6d-4b4a-8697-9926111686fc
# ╠═72e69f4a-9835-4f6b-b61d-0e5ff12319dd
# ╠═aaa31f40-a910-4f93-9bd9-601a8814a449
# ╠═8227fd27-33de-44e7-82d2-5065e5d368fb
# ╠═a08ade9f-134f-43ed-8b32-17f7750c87b2
# ╠═5d793a94-7f96-4511-841a-de622dc50662
