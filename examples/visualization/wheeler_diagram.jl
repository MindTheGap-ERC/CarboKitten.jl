# ~/~ begin <<docs/src/visualization.md#examples/visualization/wheeler_diagram.jl>>[init]
#| creates: docs/src/_fig/wheeler_diagram.png
#| requires: data/alcaps2.h5
#| collect: figures

module Script

function main()
  header, data = read_slice("data/alcaps.h5", :, 25)
  fig = wheeler_diagram(header, data)
  save("docs/src/_fig/wheeler_diagram.png", fig)
end

end

Script.main()
# ~/~ end