### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ a5f0eccc-b53d-4d41-88cd-458ec07545ab
using Pkg; Pkg.activate("../workenv")

# ╔═╡ 11b1526c-a419-445b-8b74-f5b15fba5550
using Unitful

# ╔═╡ 18d7285e-3bc5-44c1-8ba6-8dc1becefa3d
using MacroTools

# ╔═╡ f44c7bd5-a233-4605-bc1c-fe2a085cdbef
using CarboKitten.Config: Box, axes

# ╔═╡ d8b3f8a2-9110-4be0-8b86-a1d67610efe4
begin
	const Amount = typeof(1.0u"m")
	const Time = typeof(1.0u"Myr")
	const Height = typeof(1.0u"m")
	const Location = typeof(1.0u"km")
	const Rate = typeof(1.0u"m/Myr")
	const Intensity = typeof(1.0u"W/m^2")
	
	abstract type Input end
	abstract type State end
	abstract type Facies end
end

# ╔═╡ a644565f-fa1e-4040-817d-90ab36eacec9
macro spec(block)
	return QuoteNode(block)
end

# ╔═╡ bad996f3-6a72-497a-8437-b2e2259f284c
const UniformProduction = @spec begin
  using WaterDepth

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

# ╔═╡ 3777e60a-3e27-11ef-2e27-bf4856556e6d
const WaterDepth = @spec begin
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
end

# ╔═╡ 6cbc77d6-a4c9-470f-8536-a35637660a0e
struct Field
	sym::Symbol
	typ::Type
end

# ╔═╡ 374ce450-dc6d-4b4a-8697-9926111686fc
struct Struct
	mut::Bool
	fields::Vector{Field}
end

# ╔═╡ 72e69f4a-9835-4f6b-b61d-0e5ff12319dd
macro compose(cs...)
	spec = IdDict()

	pass(e::LineNumberNode) = nothing
	
	function pass(e::Expr)
		println(e)
	end
	
	for c in cs
		subspec = eval(c)
		for a in subspec.args
			pass(a)
		end
	end
	
	:()
end

# ╔═╡ aaa31f40-a910-4f93-9bd9-601a8814a449
UniformProduction.args[2].head

# ╔═╡ a08ade9f-134f-43ed-8b32-17f7750c87b2
@compose UniformProduction

# ╔═╡ Cell order:
# ╠═a5f0eccc-b53d-4d41-88cd-458ec07545ab
# ╠═11b1526c-a419-445b-8b74-f5b15fba5550
# ╠═18d7285e-3bc5-44c1-8ba6-8dc1becefa3d
# ╠═f44c7bd5-a233-4605-bc1c-fe2a085cdbef
# ╠═d8b3f8a2-9110-4be0-8b86-a1d67610efe4
# ╠═a644565f-fa1e-4040-817d-90ab36eacec9
# ╠═3777e60a-3e27-11ef-2e27-bf4856556e6d
# ╠═bad996f3-6a72-497a-8437-b2e2259f284c
# ╠═6cbc77d6-a4c9-470f-8536-a35637660a0e
# ╠═374ce450-dc6d-4b4a-8697-9926111686fc
# ╠═72e69f4a-9835-4f6b-b61d-0e5ff12319dd
# ╠═aaa31f40-a910-4f93-9bd9-601a8814a449
# ╠═a08ade9f-134f-43ed-8b32-17f7750c87b2
