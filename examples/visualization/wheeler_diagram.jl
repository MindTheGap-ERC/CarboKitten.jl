# ~/~ begin <<docs/src/visualization.md#examples/visualization/wheeler_diagram.jl>>[init]
#| creates: docs/src/_fig/wheeler_diagram.png
#| requires: data/output/alcap2.h5
#| collect: figures

module Script

using CairoMakie
using CarboKitten.Export: read_slice
using CarboKitten.Visualization: wheeler_diagram

function main()
  header, data = read_slice("data/output/alcap2.h5", :, 25)
  fig = wheeler_diagram(header, data)
  save("docs/src/_fig/wheeler_diagram.png", fig)
end

end

Script.main()
# ~/~ end
