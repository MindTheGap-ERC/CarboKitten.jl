# ~/~ begin <<docs/src/model-alcap.md#examples/model/alcap/plot.jl>>[init]

using GLMakie
using CarboKitten.Visualization

GLMakie.activate!()

save("docs/src/_fig/alcaps-alternative.png", summary_plot("data/output/alcap-example.h5"))
# ~/~ end
