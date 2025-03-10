# ~/~ begin <<docs/src/denudation/empirical.md#examples/denudation/empirical-test.jl>>[init]
#| requires: examples/denudation/empirical-test.jl
#| creates: docs/src/_fig/EmpiricalPrecipitation.png
#| collect: figures

module EmpiricalSpec
using GLMakie

using CarboKitten
using Unitful

using CarboKitten.Denudation.EmpiricalDenudationMod: empirical_denudation, slope_kernel

const slope = 30

function main()
    precip = collect(0.4:0.01:2.0)

    result = Vector{typeof(1.0u"m/Myr")}(undef,size(precip))

    for i in eachindex(precip)
        result[i] = empirical_denudation(precip[i],slope)
    end

    fig = Figure()
    ax = Axis(fig[1,1],xlabel="Precipitation (m/yr)", ylabel="Denudation rates (m/Myr)")
    lines!(ax,precip,result)
    save("docs/src/_fig/EmpiricalPrecipitation.png",fig)
end

end

EmpiricalSpec.main()
# ~/~ end
