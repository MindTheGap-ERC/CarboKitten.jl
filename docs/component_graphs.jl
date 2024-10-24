"""
The Documenter hook part was heavily inspired on the DocumenterMermaid plugin.
"""
module ComponentGraphs

using MarkdownAST: MarkdownAST, Node
using Documenter: Documenter, Selectors, DOM
using Documenter.HTMLWriter: DCtx

using GraphvizDotLang: digraph, node, edge, HTML
using MacroTools: @capture

struct ZipDefault{T,D}
	its::T
	default::D
end

function zip_default_helper(z::ZipDefault{T,D}, next::S) where {T, D, S}
	all(isnothing, next) && return nothing
	new_states = [(isnothing(i) ? nothing : i[2]) for i in next]
	results = [(isnothing(i) ? z.default : i[1]) for i in next]
	return results, new_states
end

Base.iterate(z::ZipDefault{T,D}) where {T, D} =
    zip_default_helper(z, (iterate(i) for i in z.its))

Base.iterate(z::ZipDefault{T,D}, states::S) where {T, D, S} =
    zip_default_helper(z, (isnothing(s) ? nothing : iterate(i, s)
        for (i, s) in zip(z.its, states)))

Base.eltype(z::ZipDefault{T,D}) where {T,D} = Any

Base.IteratorSize(::Type{ZipDefault{T,D}}) where {T, D} = Base.SizeUnknown()

"""
    zip_defaults(iterators...; default=nothing)

Zips its argument iterators, just like `zip` would. However, the length of the
zipped iterator is determined by the longest argument iterator. Others are
padded with the `default` value.
"""
zip_default(its...; default=nothing) = ZipDefault(its, default)

"""
    get_field_name(expr)

Takes an expression of the form `name::Type = value`, `name::Type` or `name`.
Returns the name.
"""
function get_field_name(expr)
	@capture(expr, fname_::ftype_ = fdefault_) && return fname
	@capture(expr, fname_::ftype_) && return fname
	@capture(expr, fname_) && return fname
end

format_table_row(a, b, c) = """
<tr><td align="left" border="1" sides="R"><font face="monospace" point-size="10">
$(get_field_name(a))
</font></td>
<td align="left" border="1" sides="R"><font face="monospace" point-size="10" >
$(get_field_name(b))
</font></td>
<td align="left"><font face="monospace" point-size="10">
$(get_field_name(c))
</font></td></tr>"""

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

function format_label(name, fields)
	structs = [:Input, :Facies, :State]
	all(f->f ∉ keys(fields), structs) && return "<b>$(name)</b>"
	rows = splat(format_table_row).(zip_default(
		(get(fields, s, []) for s in structs)..., default=" "))
	format_table(name, rows)
end

nodes(t, n) = union([n], (nodes(t, m) for m in t[n])...)
nodes(t) = union([], (nodes(t, n) for n in keys(t))...)
transient_dependencies(t, n) = isempty(t[n]) ?
	Set{Symbol}() :
	foldl(union, transient_dependencies.((t,), t[n]), init=t[n])

function dependency_graph(m)
	t = m.MIXIN_TREE
	g = digraph(rankdir="TD")
	for n in nodes(t)
		g |> node(string(n); shape="rect", style="rounded", margin="0.1,0.1", label=HTML(format_label(n, getfield(m, n).FIELDS)))
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

abstract type DAGExpander <: Documenter.Expanders.ExpanderPipeline end

Selectors.matcher(::Type{DAGExpander}, node, page, doc) = Documenter.iscode(node, "component-dag")
Selectors.order(::Type{DAGExpander}) = 7.8
function Selectors.runner(::Type{DAGExpander}, node, page, doc)
	# mod = Documenter.get_sandbox_module!(page.globals.meta, "component-dag", node.element.code)
	# Core.eval(mod, )
	expr = :(begin
		# using CarboKitten
		$(Meta.parse("using $(node.element.code)"))
		$(Meta.parse(split(node.element.code, ".")[end]))
	end)
	mod = eval(expr)
	io = IOBuffer()
	show(io, MIME("image/svg"), dependency_graph(mod))
	seekstart(io)
	svg = read(io, String)
	node.element = Documenter.RawNode(:html, "<div class=\"component-dag\" style=\"overflow: scroll\">$(svg)</div>")
end

end
