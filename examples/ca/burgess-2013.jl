# ~/~ begin <<docs/src/components/cellular-automata.md#examples/ca/burgess-2013.jl>>[init]
module Script
using CarboKitten
using CarboKitten.Components: CellularAutomaton as CA
using CairoMakie

function main()
  input = CA.Input(
      box = CarboKitten.Box{Periodic{2}}(
        grid_size=(50, 50), phys_scale=1.0u"m"),
      facies = fill(CA.Facies(), 3)
  )

  state = CA.initial_state(input)
  step! = CA.step!(input)

  for _ in 1:1000
    step!(state)
  end

  fig = Figure(size=(1000, 500))
  axes_indices = Iterators.flatten(eachrow(CartesianIndices((2, 4))))
  xaxis, yaxis = box_axes(input.box)
  i = 1000
  for row in 1:2
    for col in 1:4
      ax = Axis(fig[row, col], aspect=AxisAspect(1), title="step $(i)")

      if row == 2
        ax.xlabel = "x [m]"
      end
      if col == 1
        ax.ylabel = "y [m]"
      end

      heatmap!(ax, xaxis/u"m", yaxis/u"m", state.ca)
      step!(state)
      i += 1
    end
    for _ in 1:996
      step!(state)
      i += 1
    end
  end
  save("docs/src/_fig/ca-long-term.svg", fig)
end

end

Script.main()
# ~/~ end
