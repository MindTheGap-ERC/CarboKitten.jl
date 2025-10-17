# ~/~ begin <<docs/src/visualization.md#examples/visualization/wheeler_diagram.jl>>[init]

module Script

using CairoMakie
using CarboKitten.Export: read_slice
using CarboKitten.Visualization: wheeler_diagram

function main()
  header, data = read_slice("data/output/alcap-example.h5", :profile)
  fig = wheeler_diagram(header, data)
  save("docs/src/_fig/wheeler_diagram.png", fig)
end

end

Script.main()
# ~/~ end
