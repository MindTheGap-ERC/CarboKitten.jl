# ~/~ begin <<docs/src/input-methods.md#examples/tabular-sea-level/ornstein-uhlenbeck.jl>>[init]
#| creates: docs/src/_fig/OU.png
#| collect: figures

module Script

using CarboKitten

using Unitful
using CarboKitten.Components
using GLMakie
using Random

GLMakie.activate!()

const TIME_PROPERTIES = TimeProperties(
	Δt = 500u"yr",
	steps = 2000
)

function generate_ar1(mean, n, drift, variance)
"""
Arguments
- `mean`: Mean of the process
- `length`: Length of the process (number of steps)
- `drift`: Drift parameter
- `variance`: Variance of the process 
"""
    ar1 = Vector{Float64}(undef, n)
    ar1[1] = mean  # start with the mean for simplicity

    for i in 2:n
        ar1[i] = drift * ar1[i-1] + (1 - drift) * mean + randn() * sqrt(variance)
    end

    return ar1
end

const θ = 0.4 # drift
const μ = 2.0 # mean
const σ = 20 # variance

function main()

OU = generate_ar1(μ, length(time_axis(TIME_PROPERTIES)), θ, σ)

fig, ax = lines(time_axis(TIME_PROPERTIES) |> in_units_of(u"Myr"), collect(OU) .* u"m")
    ax.xlabel = "time [Myr]"
	ax.ylabel = "sea level [m]"
    save("docs/src/_fig/OU.png",fig)

end

end

Script.main()
# ~/~ end
