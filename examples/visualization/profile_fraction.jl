# ~/~ begin <<docs/src/visualization.md#examples/visualization/profile_fraction.jl>>[init]
#| creates: docs/src/_fig/profile_fraction.png
#| requires: data/output/alcap-example.h5
#| collect: figures

module Script
using GLMakie
using CarboKitten.Export: read_slice
using CarboKitten.Visualization: profile_plot!

function main()
    (header, slice) = read_slice("data/output/alcap-example.h5", :profile)
	fig = Figure()
	ax = Axis(fig[1, 1])

	x = header.axes.x
	t = header.axes.t
	
    plot = profile_plot!(x -> x[2]/sum(x), ax, header, slice; colorrange=(0, 1))
    Colorbar(fig[1, 2], plot; label=L"f_2 / f_{total}")

	save("docs/src/_fig/profile_fraction.png", fig)
	fig
end
end

Script.main()
# ~/~ end
