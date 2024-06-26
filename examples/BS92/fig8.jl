# ~/~ begin <<docs/src/bosscher-1992.md#examples/BS92/fig8.jl>>[init]
#| creates: docs/src/_fig/bs92-fig8.svg
#| requires: data/bs92-sealevel-curve.csv examples/BS92/BS92.jl
#| collect: figures

module Script
     include("BS92.jl")
     using CairoMakie

     function main()
          h0 = LinRange(0, 200, 101)
          result = hcat([BS92.model(BS92.SCENARIO_A, h).u for h in h0]...)
          t = LinRange(0, 80_000, 81)

          fig = Figure(resolution=(600,900))
          ax = Axis(fig[1,1], xlabel="initial depth (m)", ylabel="depth (m)", yreversed=true)
          for l in eachrow(result)
               lines!(ax, h0, vec(l); color=:steelblue4, linewidth=0.5)
          end
          ax = Axis(fig[2,1], xlabel="time (years)", ylabel="depth (m)", yreversed=true)
          lines!(ax, t, BS92.SCENARIO_A.sealevel(t); color=:steelblue4)

          save("docs/src/_fig/bs92-fig8.svg", fig)
     end
end

Script.main()
# ~/~ end