# ~/~ begin <<docs/src/bosscher-1992.md#src/BS92/fig8.jl>>[init]
using MindTheGap.BS92
using Plots

function main()
     h0 = LinRange(0, 200, 101)
     result = hcat([BS92.model(BS92.SCENARIO_A, h).u for h in h0]...)
     t = LinRange(0, 80_000, 81)

     plotlyjs()

     plot(h0, result',
          xaxis=("initial depth (m)"),
          yaxis=("depth (m)", :flip),
          legend_position=:none, lc=:steelblue,
          size=(700, 700), fontfamily="Merriweather,serif")

     plot!(t, BS92.SCENARIO_A.sealevel(t),
          title="sea level curve", titlelocation=:left,
          titlefontsize=12,
          xaxis=("time (years)"),
          yaxis=("depth (m)", :flip),
          guidefontsize=10,
          legend_position=:none,
          lc=:steelblue,
          inset=(1, bbox(0.11, 0.60, 0.45, 0.28)),
          subplot=2,
          framestyle=:box)

     savefig("docs/fig/bs92-fig8.html")
end
# ~/~ end