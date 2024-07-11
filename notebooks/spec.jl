### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ cbbda4e4-3f6c-11ef-05f7-6d0af577f472
macro spec(name, body)
	:(const $(esc(name)) = $(QuoteNode(body)))
end

# ╔═╡ d7096d8d-a2e8-4636-a8e0-742569eaa0e1
@spec Blah begin
	struct Hello
		a::Int
	end
end

# ╔═╡ ff0930d0-2e43-4ef0-b58c-46812cf3ca9e
Blah

# ╔═╡ Cell order:
# ╠═cbbda4e4-3f6c-11ef-05f7-6d0af577f472
# ╠═d7096d8d-a2e8-4636-a8e0-742569eaa0e1
# ╠═ff0930d0-2e43-4ef0-b58c-46812cf3ca9e
