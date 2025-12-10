# ~/~ begin <<docs/src/algorithms/stratigraphic_column.md#examples/algorithms/stratigraphic_column.jl>>[init]
module Plot

using CairoMakie
using CarboKitten.Algorithms.StratigraphicColumn: stratigraphic_column!
using Random

function moving_average(a::AbstractVector{T}, n) where {T}
	b = Vector{T}(undef, length(a))
	m = div(n, 2)
	for i in eachindex(b)
		start = min(max(i - m, 1), length(a) - n)
		b[i] = sum(a[start:start+n]) / n
	end
	return b
end

function main()
    Random.seed!(4)
    x = moving_average(randn(1000), 10) .+ 0.1

    fig = Figure()
    ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[2, 1])
    lines!(ax1, x)
    lines!(ax2, cumsum(x))

    stratigraphic_column!(x)
    lines!(ax1, x)
    lines!(ax2, cumsum(x))

    save("docs/src/_fig/monotonic_adm.svg", fig)
end

end

Plot.main()
# ~/~ end
