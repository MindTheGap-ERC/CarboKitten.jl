# ~/~ begin <<docs/src/input-methods.md#examples/tabular-sea-level/loess.jl>>[init]

module Script

using CarboKitten

using Unitful
using CarboKitten.Components

using CarboKitten.DataSets: miller_2020
using Smoothers
using GLMakie

GLMakie.activate!()

function main()

    miller_df = miller_2020()
    sort!(miller_df, [:time])

    sl = miller_df.sealevel / u"m" .|> NoUnits
    ti = miller_df.time / u"kyr" .|> NoUnits

    fig = Figure()
    ax = Axis(fig[1,1], xlabel="Time [kyr]", ylabel="Sea level [m]")
    scatter!(ax, ti, sl)
    lines!(ax, ti, Smoothers.loess(ti, sl, q = 1000)(ti); color = :tomato)
    save("docs/src/_fig/loess.png",fig)
end

end

Script.main()
# ~/~ end
