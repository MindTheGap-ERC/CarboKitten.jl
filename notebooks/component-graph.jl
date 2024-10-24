### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ df819cda-91e5-11ef-3664-4b441f84e147
using Pkg; Pkg.activate("../workenv")

# ╔═╡ d5b9a415-89fe-458f-9741-91cb245584f7
using Revise

# ╔═╡ 14c04a84-f0fb-48dd-a7f8-fd20b118ce6d
using GraphvizDotLang: digraph, node, edge, HTML

# ╔═╡ 20cab817-02e3-4658-a344-6ac8f79bcc77
using CarboKitten.Components: WaterDepth

# ╔═╡ fd12ee16-146a-4c6b-8a44-033067a65c00
using CarboKitten.Model: ALCAP2, BS92

# ╔═╡ dd686918-b698-4923-8f1e-b9b500142450
using MacroTools: @capture

# ╔═╡ 446fc567-a46e-4e66-aba5-106f189dcea8
include("../docs/component_graphs.jl")

# ╔═╡ 15e12c85-f4c4-4cb3-ab42-f3412bc967ee
g = digraph() |> edge("a", "b", "c")

# ╔═╡ c6d713d7-5c0b-49be-9c13-81f438ca35dd
ComponentGraphs.dependency_graph(ALCAP2)

# ╔═╡ 52be2bfb-fa87-41a9-8854-5fcbeecdd4b5


# ╔═╡ 890a8621-90c6-4baa-ae2b-a10454e8c140
WaterDepth.MIXIN_TREE

# ╔═╡ 023dcd9d-5657-475a-8f38-557981039d08
transient_dependencies(t, n) = isempty(t[n]) ?
	Set{Symbol}() :
	foldl(∪, transient_dependencies.((t,), t[n]), init=t[n])

# ╔═╡ b72a9f5c-1e23-4e64-9204-52674ed57ccc
transient_dependencies(ALCAP2.MIXIN_TREE, :CAProduction)

# ╔═╡ f7d0044d-8ec3-45e7-9352-00fc6cf6411f
transient_dependencies(ALCAP2.MIXIN_TREE, :Tag)

# ╔═╡ 0e72eaa4-2d92-4322-be21-ef42ffbf6b4d
nodes(t, n) = union([n], (nodes(t, m) for m in t[n])...)

# ╔═╡ b0ca2d5e-6821-4425-a3c2-4ece321c83ce
nodes(t) = union([], (nodes(t, n) for n in keys(t))...)

# ╔═╡ 4ae03c62-6dce-4e8a-a304-6ae6e44bac95
[n => getfield(ALCAP2, n).FIELDS for n in nodes(ALCAP2.MIXIN_TREE)]

# ╔═╡ cbb87393-c13c-4049-b073-7c285a7fce70
struct ZipDefault{T,D}
	its::T
	default::D
end

# ╔═╡ e332d3a5-977d-43f4-9fd6-d240080fd0ba
all(isnothing, [nothing, 1])

# ╔═╡ 174b2b09-836c-465b-9b8b-98e94a666c74
function Base.iterate(z::ZipDefault{T,D}) where {T, D}
	next = iterate.(z.its)
	all(isnothing, next) && return nothing
	new_states = [(isnothing(i) ? nothing : i[2]) for i in next]
	results = [(isnothing(i) ? z.default : i[1]) for i in next]
	return results, new_states
end

# ╔═╡ 87da86c5-6875-4b47-ac6d-cce0cdde164e
function Base.iterate(z::ZipDefault{T,D}, states::S) where {T, D, S}
	adv(i, s) = isnothing(s) ? nothing : iterate(i, s)
	next = adv.(z.its, states)
	all(isnothing, next) && return nothing
	new_states = [(isnothing(i) ? nothing : i[2]) for i in next]
	results = [(isnothing(i) ? z.default : i[1]) for i in next]
	return results, new_states
end

# ╔═╡ 56428c54-a5df-4749-a7c9-c7f04796a24f
Base.eltype(z::ZipDefault{T,D}) where {T,D} = Any

# ╔═╡ 2855e2a3-e310-4c79-aeec-7c17c9292a15
Base.IteratorSize(::Type{ZipDefault{T,D}}) where {T, D} = Base.SizeUnknown()

# ╔═╡ d800165e-2a4c-45c7-9335-d1022841dc2b
zip_default(its...; default=nothing) = ZipDefault(its, default)

# ╔═╡ 08429408-6772-48bf-8c7b-da06b92afbdc
zip_default([1], [2, 3], default=0) |> collect

# ╔═╡ f8076723-eb53-4de2-af02-d66b2727e45d
zip([1], [2, 3]) |> collect

# ╔═╡ b5e08186-08ee-44b8-a2f4-2627b7ebb8f5
get(ALCAP2.Boxes.FIELDS, :State, nothing)

# ╔═╡ 5e49f617-7529-4d97-9a0c-ed03bb3d7db0
function get_field_name(expr)
	@capture(expr, fname_::ftype_ = fdefault_) && return fname
	@capture(expr, fname_::ftype_) && return fname
	@capture(expr, fname_) && return fname
end

# ╔═╡ 7d9e95e7-425c-4559-b1aa-d56573264af4
table_row(a, b, c) = """
<tr><td align="left" border="1" sides="R"><font face="monospace" point-size="10">
$(get_field_name(a))
</font></td>
<td align="left" border="1" sides="R"><font face="monospace" point-size="10" >
$(get_field_name(b))
</font></td>
<td align="left"><font face="monospace" point-size="10">
$(get_field_name(c))
</font></td></tr>"""

# ╔═╡ 65756fee-c451-407b-8632-8083459a2e59
format_table(n, rows) = """
<table border="0" title="$(n)" cellspacing="0" cellpadding="3">
<tr><td colspan="3" border="1" sides="B"><b>$(n)</b></td></tr>
<tr><td align="left" border="1" sides="RBT" bgcolor="gray80">Input</td>
    <td align="left" border="1" sides="LRBT" bgcolor="gray80">Facies</td>
	<td align="left" border="1" sides="LBT" bgcolor="gray80">State</td>
</tr>
$(join(rows,""))
</table>
"""

# ╔═╡ 1c07f758-d042-4210-80ee-1a21397be622
function record_label(n, fields)
	structs = [:Input, :Facies, :State]
	all(f->f ∉ keys(fields), structs) && return "<b>$(n)</b>"
	rows = splat(table_row).(zip_default(
		(get(fields, s, []) for s in structs)..., default=" "))
	format_table(n, rows)
end

# ╔═╡ ae751874-3d96-409e-839c-a86b07764857
function tree_graph(m)
	t = m.MIXIN_TREE
	g = digraph(rankdir="LR")
	for n in nodes(t)
		g |> node(string(n); shape="rect", style="rounded", margin="0.1,0.1", label=HTML(record_label(n, getfield(m, n).FIELDS)))
	end
	for (k, v) in pairs(t)
		for n in v
			td = union(Set{Symbol}(), (transient_dependencies(t, m) for m in setdiff(v, [n]))...)
			if n ∉ td
				g |> edge(string(n), string(k))
			end
		end
	end
	g
end

# ╔═╡ ed64dc59-8143-4e89-821c-a94a532ad97c
tree_graph(WaterDepth)

# ╔═╡ cf688a1b-b487-4ef3-ba40-5e202ad40291
tree_graph(ALCAP2)

# ╔═╡ 571be0b7-3dee-48b7-963d-45688d72280e
tree_graph(ALCAP2.ActiveLayer)

# ╔═╡ aca3c5c8-8e82-4e9d-a2f0-b45869dd59d1
tree_graph(BS92)

# ╔═╡ Cell order:
# ╠═df819cda-91e5-11ef-3664-4b441f84e147
# ╠═d5b9a415-89fe-458f-9741-91cb245584f7
# ╠═14c04a84-f0fb-48dd-a7f8-fd20b118ce6d
# ╠═15e12c85-f4c4-4cb3-ab42-f3412bc967ee
# ╠═446fc567-a46e-4e66-aba5-106f189dcea8
# ╠═c6d713d7-5c0b-49be-9c13-81f438ca35dd
# ╠═20cab817-02e3-4658-a344-6ac8f79bcc77
# ╠═fd12ee16-146a-4c6b-8a44-033067a65c00
# ╠═52be2bfb-fa87-41a9-8854-5fcbeecdd4b5
# ╠═890a8621-90c6-4baa-ae2b-a10454e8c140
# ╠═023dcd9d-5657-475a-8f38-557981039d08
# ╠═b72a9f5c-1e23-4e64-9204-52674ed57ccc
# ╠═f7d0044d-8ec3-45e7-9352-00fc6cf6411f
# ╠═0e72eaa4-2d92-4322-be21-ef42ffbf6b4d
# ╠═b0ca2d5e-6821-4425-a3c2-4ece321c83ce
# ╠═4ae03c62-6dce-4e8a-a304-6ae6e44bac95
# ╠═cbb87393-c13c-4049-b073-7c285a7fce70
# ╠═e332d3a5-977d-43f4-9fd6-d240080fd0ba
# ╠═174b2b09-836c-465b-9b8b-98e94a666c74
# ╠═87da86c5-6875-4b47-ac6d-cce0cdde164e
# ╠═56428c54-a5df-4749-a7c9-c7f04796a24f
# ╠═2855e2a3-e310-4c79-aeec-7c17c9292a15
# ╠═d800165e-2a4c-45c7-9335-d1022841dc2b
# ╠═08429408-6772-48bf-8c7b-da06b92afbdc
# ╠═f8076723-eb53-4de2-af02-d66b2727e45d
# ╠═b5e08186-08ee-44b8-a2f4-2627b7ebb8f5
# ╠═dd686918-b698-4923-8f1e-b9b500142450
# ╠═5e49f617-7529-4d97-9a0c-ed03bb3d7db0
# ╠═7d9e95e7-425c-4559-b1aa-d56573264af4
# ╠═65756fee-c451-407b-8632-8083459a2e59
# ╠═1c07f758-d042-4210-80ee-1a21397be622
# ╠═ae751874-3d96-409e-839c-a86b07764857
# ╠═ed64dc59-8143-4e89-821c-a94a532ad97c
# ╠═cf688a1b-b487-4ef3-ba40-5e202ad40291
# ╠═571be0b7-3dee-48b7-963d-45688d72280e
# ╠═aca3c5c8-8e82-4e9d-a2f0-b45869dd59d1
