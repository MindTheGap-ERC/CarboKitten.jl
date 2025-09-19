# ~/~ begin <<docs/src/denudation/chemical.md#examples/denudation/dissolution-test.jl>>[init]
#| requires: examples/denudation/dissolution-test.jl
#| creates:
#|   - docs/src/_fig/KHTemp.png
#|   - docs/src/_fig/Equilibrium_Concs.png
#|   - docs/src/_fig/DissolutionExample.png
#| collect: figures

module DissolutionSpec
using CairoMakie
using Unitful
using CarboKitten.Denudation.DissolutionMod: equilibrium, dissolution, karst_denudation_parameters
using CarboKitten: Box
using CarboKitten.Stencil: Periodic, Reflected, stencil


@kwdef struct facies
    infiltration_coefficient :: Float64
    mass_density :: typeof(1.0u"kg/m^3")
    reactive_surface :: typeof(1.0u"m^2/m^3")
end


@kwdef struct Dissolution
    temp :: Any
    precip :: Float64
    pco2 :: Float64
    reactionrate :: Float64
end

@kwdef struct state
    ca::Matrix{Int}
end

const BOX = Box{Periodic{2}}(grid_size=(5, 5), phys_scale=1.0u"km")


const Facies = facies(infiltration_coefficient = 0.5,
mass_density = 2.73u"kg/m^3",
reactive_surface = 1000u"m^2/m^3")

const DIS = Dissolution(temp = 285,
precip  = 1.0,
pco2 = 10^(-1.5),
reactionrate = 0.1)

const STATE = state(ca = 
[ 0  0  1  0  0
0  1  0  1  1
0  0  1  0  1
1  0  1  1  0
1  1  1  1  0])

const WD = -10 .* [ 0.0663001  0.115606  0.646196
0.601523   0.130196  0.390821
0.864734   0.902935  0.670354]

function main()
    temp = collect(293:0.5:303)

    Eq_result = Array{Any,1}(undef,length(temp))
    Para_result = Array{Any,1}(undef,length(temp))
    Dis_result = Array{Any,1}(undef,length(temp))
    Dis_result_wd = Array{Any,1}(undef,length(WD))

    for i in eachindex(temp)
        Para_result[i] = karst_denudation_parameters(temp[i])
    end
    Para_result_KH = [x.KH for x in Para_result]

    for i in eachindex(temp)
       Eq_result[i] =  equilibrium(temp[i],DIS.pco2,DIS.precip,Facies)
    end
    Eq_result_concs = [x.concentration for x in Eq_result]
    for i in eachindex(temp)
        Dis_result[i] = dissolution(temp[i],DIS.precip,DIS.pco2,DIS.reactionrate,WD[1],Facies)
    end
    Dis_result
    for i in eachindex(WD)
        Dis_result_wd[i] = dissolution(temp[1],DIS.precip,DIS.pco2,DIS.reactionrate,WD[i],Facies)
    end
    Dis_result_wd

    fig1 = Figure()
    ax1 = Axis(fig1[1,1],xlabel="Temp (K)", ylabel=" concentration (mol/L)")
    lines!(ax1,temp,Eq_result_concs)
    save("docs/src/_fig/Equilibrium_Concs.png",fig1)

    fig2 = Figure()
    ax2 = Axis(fig2[1,1],xlabel="Temp (K)", ylabel="KH")
    lines!(ax2,temp,Para_result_KH)
    save("docs/src/_fig/KHTemp.png",fig2)

    fig3 = Figure()
    ax3 = Axis(fig3[1,1],xlabel="Temp (K)", ylabel="Denudation rates (m/Myr)")
    lines!(ax3,temp,Dis_result)

    ax4 = Axis(fig3[1,2],xlabel="Water Depth (m)", ylabel="Denudation rates (m/Myr)")
    scatter!(ax4,vec(WD),Dis_result_wd)
    fig3
    save("docs/src/_fig/DissolutionExample.png",fig3)
end
end

DissolutionSpec.main()
# ~/~ end