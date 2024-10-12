### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 8433507c-88b1-11ef-201c-df0d6533a9fe
using Pkg; Pkg.activate("../workenv")

# ╔═╡ af8f7212-a873-4d76-8a7a-e876c1ebb6e5
using Revise

# ╔═╡ e671c760-0c57-432b-8cc1-4c7e85b8b514
using GLMakie

# ╔═╡ e3df401a-5885-4a77-812a-aa739ca557b6
using CarboKitten.Export: read_slice, Header, DataSlice

# ╔═╡ 9ddaa447-5f7c-4293-b5bc-668857d8ad16
using CarboKitten.Visualization

# ╔═╡ 3246dc1b-bfd1-4f37-8f49-3e467907c1ac
using CarboKitten.Utility: in_units_of

# ╔═╡ 9e7add76-7d16-4bf0-b9d0-4a81b28865cd
using Unitful

# ╔═╡ cee9e8a5-d897-42ba-a2b1-53bee950e366
let
	header, data = read_slice("../data/output/alcap1.h5", :, 25)
	n_facies = size(data.production)[1]
	fig = Figure(size=(1000, 1000))
	ax1 = Axis(fig[1,1:2])
	sediment_profile!(ax1, header, data)
	ax2 = Axis(fig[3,1])
	ax3 = Axis(fig[3,2])
	sm, df = wheeler_diagram!(ax2, ax3, header, data)
	Colorbar(fig[2,1], sm; vertical=false, label="sediment accumulation [m/Myr]")
	Colorbar(fig[2,2], df; vertical=false, label="dominant facies", ticks=1:n_facies)

	ax4 = Axis(fig[3,3], title="sealevel curve", ylabel="time [Myr]", xlabel="sealevel [m]", limits=(nothing, (0.0, header.axes.t[end] |> in_units_of(u"Myr"))))
	lines!(ax4, header.sea_level |> in_units_of(u"m"), header.axes.t |> in_units_of(u"Myr"))
	fig
end

# ╔═╡ 36d424f2-dd37-4fbc-84c4-a0c664ce83b2
fieldnames(Header)

# ╔═╡ Cell order:
# ╠═8433507c-88b1-11ef-201c-df0d6533a9fe
# ╠═af8f7212-a873-4d76-8a7a-e876c1ebb6e5
# ╠═e671c760-0c57-432b-8cc1-4c7e85b8b514
# ╠═e3df401a-5885-4a77-812a-aa739ca557b6
# ╠═9ddaa447-5f7c-4293-b5bc-668857d8ad16
# ╠═3246dc1b-bfd1-4f37-8f49-3e467907c1ac
# ╠═9e7add76-7d16-4bf0-b9d0-4a81b28865cd
# ╠═cee9e8a5-d897-42ba-a2b1-53bee950e366
# ╠═36d424f2-dd37-4fbc-84c4-a0c664ce83b2
